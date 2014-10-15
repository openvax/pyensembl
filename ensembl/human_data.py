from glob import glob
import logging
from os.path import join, exists, split
from os import remove

from gtf import load_gtf_as_dataframe
from locus import normalize_chromosome

import datacache
import numpy as np
import pandas as pd

MIN_ENSEMBL_RELEASE = 48
MAX_ENSEMBL_RELEASE = 77



def _check_release(release):
    """
    Convert a user-provided release number into
    an integer, check to make sure it's in the
    valid range of Ensembl releases
    """
    try:
        release = int(release)
    except:
       assert False, "%s is not a valid Ensembl release" % release
    assert release >= MIN_ENSEMBL_RELEASE
    assert release <= MAX_ENSEMBL_RELEASE
    return release



# mapping from Ensembl release to which reference assembly it uses
_human_references = {}

# Ensembl release 48-54 use NCBI36 as a reference
for i in xrange(48,55):
    _human_references[i] = 'NCBI36'

# Ensembl releases 55-75 use CRCh37 as a reference
for i in xrange(55,76):
    _human_references[i] = 'GRCh37'

# Ensembl releases 76 and 77 use GRCh38
for i in xrange(76,78):
    _human_references[i] = 'GRCh38'

def _which_human_reference(release):
    release = _check_release(release)
    assert release in _human_references, \
        "No reference found for release %d" % release
    return _human_references[release]


# directory which contains GTF files, missing the release number
URL_DIR_TEMPLATE = 'ftp://ftp.ensembl.org/pub/release-%d/gtf/homo_sapiens/'
FILENAME_TEMPLATE = "Homo_sapiens.%s.%d.gtf.gz"

class EnsemblRelease(object):
    def __init__(self, release):
        self.release = _check_release(release)
        self.gtf_url_dir = URL_DIR_TEMPLATE % self.release
        self.reference_name =  _which_human_reference(self.release)
        self.gtf_filename = FILENAME_TEMPLATE  % (
            self.reference_name, self.release
        )
        self.gtf_url = join(self.gtf_url_dir, self.gtf_filename)
        self.clear_cached_objects()

    def clear_cached_objects(self):
        """
        Reset all the cached fields on this object, forcing
        everything to be reloaded from disk
        """
         # lazily download GTF data if anything beyond path/URL is necessary
        self._local_gtf_path = None

        # lazily load DataFrame if necessary
        self._df = None

        # lazily cache DataFrame with expanded attribute columns
        self._local_csv_path = None

        # if we load a subset of the entries for a specific contig
        # cache it by the contig name in this dictionary
        self._contig_dfs = {}

        # paths to partial DataFrames for particular contigs, saved as CSVs
        self._contig_csv_paths = {}


    def base_gtf_filename(self):
        """
        Trim extensions such as ".gtf" or ".gtf.gz",
        leaving only the base filename which should be used
        to construct other derived filenames for cached data.
        """
        assert ".gtf" in self.gtf_filename, \
            "GTF filename must contain .gtf extension: %s" % self.gtf_filename
        parts = self.gtf_filename.split(".gtf")
        return parts[0]


    def delete_cached_files(self):
        base = self.base_gtf_filename()
        dirpath = self.local_gtf_dir()
        for path in glob(join(dirpath, base + "*")):
            logging.info("Deleting cached file %s", path)
            remove(path)
        self.clear_cached_objects()

    def local_gtf_path(self):
        """
        Returns local path to GTF file for given release of Ensembl,
        download from the Ensembl FTP server if not already cached.
        """
        if self._local_gtf_path is None:
            self._local_gtf_path = datacache.fetch_file(
                self.gtf_url,
                filename=self.gtf_filename,
                decompress=False,
                subdir="ensembl")
        assert self._local_gtf_path
        return self._local_gtf_path

    def local_gtf_dir(self):
        return split(self.local_gtf_path())[0]

    def local_csv_path(self):
        """
        Path to CSV which the annotation data with expanded columns
        for optional attributes
        """
        if self._local_csv_path is None:
            gtf_path = self.local_gtf_path()
            base = self.base_gtf_filename()
            csv_filename = base + ".expanded.csv"
            dirpath = self.local_gtf_dir()
            self._local_csv_path = join(dirpath, csv_filename)
        return self._local_csv_path

    def local_csv_path_for_contig(self, contig_name):
        """
        Path to CSV file containing subset of Ensembl data
        restricted to given contig_name
        """
        dirpath = self.local_gtf_dir()
        base = self.base_gtf_filename()
        contig_name = normalize_chromosome(contig_name)
        return join(dirpath, base + ".contig.%s.csv" % contig_name)


    def _load_from_gtf(self):
        """
        Parse this release's GTF file and load it as a Pandas DataFrame
        """
        path = self.local_gtf_path()
        return load_gtf_as_dataframe(path)

    def _load_from_csv(self):
        """
        If we've already saved the DataFrame as a CSV
        (with the attributes field expanded into columns),
        then load it. Otherwise parse the GTF file, and save it
        as a CSV
        """
        csv_path = self.local_csv_path()
        if exists(csv_path):
            print "Reading Dataframe from %s" % csv_path
            df = pd.read_csv(csv_path)
            # by default, Pandas will infer the type as int,
            # then switch to str when it hits non-numerical
            # chromosomes. Make sure whole column has the same type
            df['seqname'] = df['seqname'].map(str)
            return df
        else:
            df = self._load_from_gtf()
            print "Saving expanded DataFrame to %s" % csv_path
            df.to_csv(csv_path, index=False)
            return df

    def dataframe(self):
        if self._df is None:
            self._df = self._load_from_csv()
        assert self._df is not None
        return self._df

    def dataframe_for_contig(self, contig_name):
        """
        Load a subset of the Ensembl data for a specific contig
        """
        contig_name = normalize_chromosome(contig_name)
        if contig_name not in self._contig_dfs:
            contig_csv_path = self.local_csv_path_for_contig(contig_name)
            df = self.dataframe()
            mask = df.seqname == contig_name
            subset = df[mask]
            assert len(subset) > 0, "Contig not found: %s" % contig_name
            self._contig_dfs[contig_name] = subset
        return self._contig_dfs[contig_name]

    def dataframe_at_loci(self, contig_name, start, stop):
        """
        Subset of entries which overlap an inclusive range of loci
        """
        contig_name = normalize_chromosome(contig_name)
        df_chr = self.dataframe_for_contig(contig_name)

        # find genes whose start/end boundaries overlap with the position
        overlap_start = df_chr.start <= stop
        overlap_end = df_chr.end >= start
        overlap = overlap_start & overlap_end
        df_overlap = df_chr[overlap]
        return df_overlap

    def dataframe_at_locus(self, contig_name, position):
        return self.dataframe_at_loci(
            contig_name=contig_name,
            start=position,
            stop=position)

    def property_values_at_loci(self, contig_name, start, stop, property_name):
        df = self.dataframe_at_loci(contig_name, start, stop)
        assert property_name in df, \
            "Unknown Ensembl property: %s" % property_name
        # drop empty strings and NaNs
        col = [
            x for  x in df[property_name]
            if x and (not np.isreal(x) or not np.isnan(x))
        ]
        return list(sorted(set(col)))

    def property_values_at_locus(self, contig_name, position, property_name):
        return self.property_values_at_loci(
            contig_name=contig_name,
            start=position,
            stop=position,
            property_name=property_name)

    def gene_ids_at_locus(self, contig_name, position):
        return self.property_values_at_locus(contig_name, position, 'gene_id')

    def gene_names_at_locus(self, contig_name, position):
        return self.property_values_at_locus(contig_name, position, 'gene_name')

    def exon_ids_at_locus(self, contig_name, position):
        return self.property_values_at_locus(contig_name, position, 'exon_id')

    def transcript_ids_at_locus(self, contig_name, position):
        return self.property_values_at_locus(
            contig_name, position, 'transcript_id')

    def transcript_names_at_locus(self, contig_name, position):
        return self.property_values_at_locus(
            contig_name, position, 'transcript_name')

    def gene_ids_at_loci(self, contig_name, start, stop):
        return self.property_values_at_loci(
            contig_name=contig_name,
            start=start,
            stop=stop,
            property_name='gene_id')

    def gene_names_at_loci(self, contig_name, start, stop):
        return self.property_values_at_loci(
            contig_name=contig_name,
            start=start,
            stop=stop,
            property_name='gene_name')

    def exon_ids_at_loci(self, contig_name, start, stop):
        return self.property_values_at_loci(
            contig_name=contig_name,
            start=start,
            stop=stop,
            property_name='exon_id')

    def transcript_ids_at_loci(self, contig_name, start, stop):
        return self.property_values_at_loci(
            contig_name=contig_name,
            start=start,
            stop=stop,
            property_name='transcript_id')

    def transcript_names_at_loci(self, contig_name, start, stop):
        return self.property_values_at_loci(
            contig_name=contig_name,
            start=start,
            stop=stop,
            property_name='transcript_name')


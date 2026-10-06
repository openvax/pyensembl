"""Explicit assembly/date selections from the new Ensembl FTP platform."""

from datetime import datetime
import re

from .ensembl_url_templates import ENSEMBL_PLATFORM_FTP_SERVER, make_dated_release_urls
from .genome import Genome
from .species import find_species_by_name


class EnsemblAnnotation(Genome):
    """A dated geneset for one assembly on the new Ensembl platform.

    Select an existing assembly accession, provider and annotation date from
    Ensembl's downloads page. These are not numbered ``EnsemblRelease`` values
    or the website's YYYY-MM release label. Construction never downloads data;
    call ``download()`` then ``index()`` to install the selected files.

    Parameters
    ----------
    assembly_accession : str
        Versioned INSDC accession, e.g. GCA_000001405.29.
    annotation_date : str
        Annotation directory in YYYY_MM form, e.g. 2023_03.
    provider : str
        Existing provider directory, usually "ensembl" or "community".
    include_alt : bool
        Select genes-including_alt.gtf.gz instead of genes.gtf.gz. Check that
        this file exists for the selected dataset; coverage is never inferred.
    genome_fasta : bool
        Also install combined reference DNA. False keeps DNA optional.
    genome_fasta_mask : str
        "none", "soft" or "hard", selecting the corresponding genome FASTA.
    species : str, optional
        Known species name, used to guard species-specific alias lookups.
        The assembly accession remains the authoritative dataset selection.
    reference_name : str, optional
        Display name for the assembly; defaults to the accession.
    cache_directory_path : str, optional
        Explicit cache directory. Use a distinct directory per dataset.
    server : str
        Root of a mirror using the same accession/provider/date layout.
    """

    def __init__(
        self, assembly_accession, annotation_date, *, provider="ensembl",
        include_alt=False, genome_fasta=False, genome_fasta_mask="none",
        species=None, reference_name=None, cache_directory_path=None,
        server=ENSEMBL_PLATFORM_FTP_SERVER,
    ):
        match = re.fullmatch(r"(GC[AF])_(\d{3})(\d{3})(\d{3})\.(\d+)", assembly_accession)
        if match is None or int(match.group(5)) < 1:
            raise ValueError("Expected a versioned GCA/GCF assembly accession")
        if not re.fullmatch(r"\d{4}_\d{2}", annotation_date):
            raise ValueError("annotation_date must use YYYY_MM")
        try:
            datetime.strptime(annotation_date, "%Y_%m")
        except ValueError as error:
            raise ValueError("annotation_date must be a valid YYYY_MM date") from error
        if not re.fullmatch(r"[a-z][a-z0-9_]*", provider):
            raise ValueError("provider must be one lowercase directory name")
        if not isinstance(include_alt, bool) or not isinstance(genome_fasta, bool):
            raise TypeError("include_alt and genome_fasta must be bool values")
        self.assembly_accession = assembly_accession
        self.annotation_date = annotation_date
        self.provider = provider
        self.include_alt = include_alt
        self._genome_fasta_option = genome_fasta
        self.genome_fasta_mask = genome_fasta_mask
        self.species = find_species_by_name(species) if species is not None else None
        self.server = server.rstrip("/")
        urls = make_dated_release_urls(
            assembly_accession, provider, annotation_date, include_alt=include_alt,
            genome_fasta_mask=genome_fasta_mask, server=self.server,
        )
        self.download_url = urls.directory
        coverage = "including-alt" if include_alt else "primary"
        super().__init__(
            reference_name=reference_name or assembly_accession,
            annotation_name="ensembl-%s-%s-%s-" % (provider, assembly_accession, coverage),
            annotation_version=annotation_date,
            gtf_path_or_url=urls.gtf,
            transcript_fasta_paths_or_urls=[urls.cdna],
            protein_fasta_paths_or_urls=[urls.pep],
            genome_fasta_path_or_url=urls.genome_fasta if genome_fasta else None,
            cache_directory_path=cache_directory_path,
        )

    def to_dict(self):
        """Preserve dataset selection and optional DNA during serialization."""
        return dict(
            assembly_accession=self.assembly_accession,
            annotation_date=self.annotation_date, provider=self.provider,
            include_alt=self.include_alt, genome_fasta=self._genome_fasta_option,
            genome_fasta_mask=self.genome_fasta_mask,
            species=self.species.latin_name if self.species is not None else None,
            reference_name=self.reference_name, cache_directory_path=self.cache_directory_path,
            server=self.server,
        )

    @classmethod
    def from_dict(cls, state_dict):
        return cls(**state_dict)

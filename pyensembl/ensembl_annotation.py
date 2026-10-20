"""Compatibility import for the original explicit assembly/date constructor."""

from .ensembl_release import EnsemblRelease
from .ensembl_url_templates import ENSEMBL_PLATFORM_FTP_SERVER
from .ensembl_versions import is_dated_release
from .genome import Genome


class EnsemblAnnotation(EnsemblRelease):
    """Compatibility wrapper; new code should use :class:`EnsemblRelease`.

    Replace ``EnsemblAnnotation(accession, date, **options)`` with
    ``EnsemblRelease(date, assembly_accession=accession, **options)``.
    This wrapper retains its original defaults, cache identity, and serialized
    class path. It emits no warning, including when loading old saved objects.
    """

    def __init__(
        self, assembly_accession, annotation_date, *, provider="ensembl",
        include_alt=False, genome_fasta=False, genome_fasta_mask="none",
        species=None, reference_name=None, cache_directory_path=None,
        server=ENSEMBL_PLATFORM_FTP_SERVER, _resolved=None,
    ):
        if not isinstance(genome_fasta, bool) or not isinstance(include_alt, bool):
            raise TypeError("include_alt and genome_fasta must be bool")
        if not isinstance(annotation_date, str) or not is_dated_release(annotation_date):
            raise ValueError("annotation_date must use a valid YYYY_MM date")
        super().__init__(
            annotation_date, species=species, server=server,
            assembly_accession=assembly_accession, provider=provider,
            include_alt=include_alt, genome_fasta=genome_fasta,
            genome_fasta_mask=genome_fasta_mask, reference_name=reference_name,
            cache_directory_path=cache_directory_path,
            _resolved=_resolved, _legacy_annotation=True,
        )

    # Preserve the original identity for saved sets/dictionaries.
    __eq__ = Genome.__eq__
    __hash__ = Genome.__hash__
    __str__ = Genome.__str__

    def to_dict(self):
        return dict(
            assembly_accession=self.assembly_accession,
            annotation_date=self.annotation_date, provider=self.provider,
            include_alt=self.include_alt, genome_fasta=bool(self._genome_fasta_option),
            genome_fasta_mask=self.genome_fasta_mask,
            species=self.species.latin_name if self.species is not None else None,
            reference_name=self.reference_name, cache_directory_path=self.cache_directory_path,
            server=self.server, _resolved=self._resolved,
        )

    @classmethod
    def from_dict(cls, state_dict):
        return cls(**state_dict)

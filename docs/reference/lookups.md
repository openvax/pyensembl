# Lookup methods

These methods apply to `EnsemblRelease`, `EnsemblAnnotation` and `Genome`.
The [complete Genome reference](genome.md) includes signatures and all public methods.

The `EnsemblRelease` object has methods to let you access all possible
combinations of the annotation features `gene_name`, `gene_id`,
`transcript_name`, `transcript_id`, `exon_id` as well as the location of
these genomic elements (contig, start position, end position, strand).

## Genes

<dl>
<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20genes(">genes(contig=None, strand=None, biotype=None)</a></dt>
<dd>Returns a list of <a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/gene.py#:~:text=class%20Gene(">Gene</a> objects, optionally restricted to a particular contig,
strand, or <code>gene_biotype</code>.</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20genes_at_locus(">genes_at_locus(contig, position, end=None, strand=None)</a></dt>
<dd>Returns a list of <a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/gene.py#:~:text=class%20Gene(">Gene</a> objects overlapping a particular position on a contig,
optionally extend into a range with the end parameter and restrict to
forward or backward strand by passing strand='+' or strand='-'.</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20gene_by_id(">gene_by_id(gene_id)</a></dt>
<dd>Return a <a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/gene.py#:~:text=class%20Gene(">Gene</a> object for given Ensembl gene ID (e.g. "ENSG00000068793").</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20gene_names(">gene_names(contig=None, strand=None)</a></dt>
<dd>Returns all gene names in the annotation database, optionally restricted
to a particular contig or strand.</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20genes_by_name(">genes_by_name(gene_name)</a></dt>
<dd>Get all the unique genes with the given name (there might be multiple
due to copies in the genome), return a list containing a <a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/gene.py#:~:text=class%20Gene(">Gene</a> object for each
distinct ID.</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20gene_by_protein_id(">gene_by_protein_id(protein_id)</a></dt>
<dd>Find <a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/gene.py#:~:text=class%20Gene(">Gene</a> associated with the given Ensembl protein ID (e.g. "ENSP00000350283")</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20gene_names_at_locus(">gene_names_at_locus(contig, position, end=None, strand=None)</a></dt>
<dd>Names of genes overlapping with the given locus, optionally restricted by strand.
(returns a list to account for overlapping genes)</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20gene_name_of_gene_id(">gene_name_of_gene_id(gene_id)</a></dt>
<dd>Returns name of gene with given gene ID.</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20gene_name_of_transcript_id(">gene_name_of_transcript_id(transcript_id)</a></dt>
<dd>Returns name of gene associated with given transcript ID.</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20gene_name_of_transcript_name(">gene_name_of_transcript_name(transcript_name)</a></dt>
<dd>Returns name of gene associated with given transcript name.</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20gene_name_of_exon_id(">gene_name_of_exon_id(exon_id)</a></dt>
<dd>Returns name of gene associated with given exon ID.</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20gene_ids(">gene_ids(contig=None, strand=None, biotype=None)</a></dt>
<dd>Return all gene IDs in the annotation database, optionally restricted by
chromosome name, strand, or <code>gene_biotype</code>.</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20gene_ids_of_gene_name(">gene_ids_of_gene_name(gene_name)</a></dt>
<dd>Returns all Ensembl gene IDs with the given name.</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20nearest_gene(">nearest_gene(contig, position, end=None, strand=None)</a></dt>
<dd>Returns <code>(distance, <a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/gene.py#:~:text=class%20Gene(">Gene</a>)</code> for the gene whose locus is nearest to the
position (or position..end interval) on the given contig — even when no
gene overlaps. Returns <code>(inf, None)</code> when no candidates exist.</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20merged_gene_intervals(">merged_gene_intervals(contig, strand=None)</a></dt>
<dd>Returns the union of all gene loci on the contig as a sorted list of
non-overlapping <code>(start, end)</code> tuples. Adjacent intervals
(<code>end+1 == next start</code>) are merged into one.</dd>

</dl>

## Transcripts

<dl>
<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20transcripts(">transcripts(contig=None, strand=None, biotype=None)</a></dt>
<dd>Returns a list of <a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/transcript.py#:~:text=class%20Transcript(">Transcript</a> objects for all transcript entries in the
Ensembl database, optionally restricted to a particular contig, strand, or
<code>transcript_biotype</code>.</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20transcript_by_id(">transcript_by_id(transcript_id)</a></dt>
<dd>Construct a <a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/transcript.py#:~:text=class%20Transcript(">Transcript</a> object for given Ensembl transcript ID (e.g. "ENST00000369985")</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20transcripts_by_name(">transcripts_by_name(transcript_name)</a></dt>
<dd>Returns a list of <a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/transcript.py#:~:text=class%20Transcript(">Transcript</a> objects for every transcript matching the given name.</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20transcript_names(">transcript_names(contig=None, strand=None)</a></dt>
<dd>Returns all transcript names in the annotation database.</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20transcript_ids(">transcript_ids(contig=None, strand=None, biotype=None)</a></dt>
<dd>Returns all transcript IDs in the annotation database.</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20transcript_ids_of_gene_id(">transcript_ids_of_gene_id(gene_id)</a></dt>
<dd>Return IDs of all transcripts associated with given gene ID.</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20transcript_ids_of_gene_name(">transcript_ids_of_gene_name(gene_name)</a></dt>
<dd>Return IDs of all transcripts associated with given gene name.</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20transcript_ids_of_transcript_name(">transcript_ids_of_transcript_name(transcript_name)</a></dt>
<dd>Find all Ensembl transcript IDs with the given name.</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20transcript_ids_of_exon_id(">transcript_ids_of_exon_id(exon_id)</a></dt>
<dd>Return IDs of all transcripts associated with given exon ID.</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20nearest_transcript(">nearest_transcript(contig, position, end=None, strand=None)</a></dt>
<dd>Returns <code>(distance, <a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/transcript.py#:~:text=class%20Transcript(">Transcript</a>)</code> to the closest transcript on the contig.
Returns <code>(inf, None)</code> when no candidates exist.</dd>
</dl>

## Exons

<dl>
<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20exon_ids(">exon_ids(contig=None, strand=None)</a></dt>
<dd>Returns a list of exon IDs in the annotation database, optionally restricted
by the given chromosome and strand.</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20exon_by_id(">exon_by_id(exon_id)</a></dt>
<dd>Construct an <a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/exon.py#:~:text=class%20Exon(">Exon</a> object for given Ensembl exon ID (e.g. "ENSE00001209410")</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20exon_ids_of_gene_id(">exon_ids_of_gene_id(gene_id)</a></dt>
<dd>Returns a list of exon IDs associated with a given gene ID.</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20exon_ids_of_gene_name(">exon_ids_of_gene_name(gene_name)</a></dt>
<dd>Returns a list of exon IDs associated with a given gene name.</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20exon_ids_of_transcript_id(">exon_ids_of_transcript_id(transcript_id)</a></dt>
<dd>Returns a list of exon IDs associated with a given transcript ID.</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20exon_ids_of_transcript_name(">exon_ids_of_transcript_name(transcript_name)</a></dt>
<dd>Returns a list of exon IDs associated with a given transcript name.</dd>
</dl>

## Reference DNA

These need reference DNA; see [Reference DNA](../guides/reference-dna.md).

<dl>
<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20sequence(">sequence(contig, start, end, mask="upper", strand="+")</a></dt>
<dd>Returns the bases from <code>start</code> to <code>end</code> (one-based, inclusive) on the plus strand, or their reverse complement with <code>strand="-"</code>. <code>mask="raw"</code> keeps soft-masked lowercase.</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20download_genome_fasta(">download_genome_fasta(overwrite=False, show_progress=False)</a></dt>
<dd>Downloads the configured reference DNA without annotation data; does nothing if it is already installed. <code>show_progress=True</code> shows progress bars.</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20index_genome_fasta(">index_genome_fasta(overwrite=False, show_progress=False)</a></dt>
<dd>Builds the DNA index now rather than on the first query.</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20fasta(">fasta</a></dt>
<dd>A <a href="https://github.com/mdshw5/pyfaidx">pyfaidx</a> reader with zero-based, half-open slices, or <code>None</code> when DNA is not configured or not installed.</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20genome_fasta_path(">genome_fasta_path</a></dt>
<dd>Path of the installed, uncompressed FASTA, or <code>None</code>.</dd>
</dl>

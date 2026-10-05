# Documentation style

Write for someone selecting reference data, running a query and interpreting
the result. Follow the [mhctools writing guide](https://github.com/openvax/mhctools/blob/fd4f3d6/docs/dev/documentation-style.md)
and its [quiet typography](https://github.com/openvax/mhctools/blob/fd4f3d6/docs/stylesheets/docs.css).
The site uses [MkDocs](https://www.mkdocs.org/) with
[Material](https://squidfunk.github.io/mkdocs-material/); API pages use
[mkdocstrings](https://mkdocstrings.github.io/).

## Structure

Put installation and the first useful examples directly on the homepage.
Readers should reach a result without having to find a separate tutorial.
Use one explicitly selected species, assembly and annotation version.
Explain IDs, coordinates, strand and sequence interpretation alongside the
output. Task guides cover selection, custom files, aliases, caches and DNA.
Reference pages give complete interfaces and metadata. Maintainer workflows
belong under Development.

Use descriptive headings and label guides by the task a reader wants to do.
Keep navigation groups collapsed until they are needed. Split pages for
substantial tasks or reference material, not for each step of a short example.
Order common lookups and sequence tasks before reference selection, cache
management and custom-data setup. Use the same reference for the main examples;
put examples requiring another species or release after the main workflow.
Avoid repeating the sidebar in the page body.
Comparison tables should have short entries; place long explanations below
them or on linked reference pages.

## Prose and literal syntax

Use ordinary text for biological concepts, units and library names. Link names
on first use in a section and in comparison tables. Reserve inline code for
literal API names, parameters, fields, commands, paths and values.

| Use code | Use ordinary text |
| --- | --- |
| `genes_by_name()`, `EnsemblRelease` | [Ensembl](https://www.ensembl.org/), annotation, genome assembly |
| `strand`, `gene_id`, `aliases=` | Chromosome, transcript, exon, base, amino acid |
| `"GRCh38"`, `93`, `None` | Human GRCh38, Ensembl release 93 |
| `pyensembl install --release 93` | Installation, indexing, cache readiness |

Keep scientific qualifications next to their claims: assembly and annotation
selection, one-based inclusive coordinates, ambiguity, source provenance and
the difference between spliced cDNA and genomic DNA. Matching file names or
contig names alone do not prove assembly identity.

## Examples and presentation

Prefer complete examples with imports and real inputs. State required downloads
before queries. Pin data and show output when it helps interpretation. Label
templates and replaceable paths; do not present them as immediately runnable.
Keep working examples, source links, page URLs and linked anchors where practical.
README compatibility anchors point to the migrated guidance.

Use system fonts, a readable line length and restrained headings. Wide tables
and code blocks may scroll within their containers. Use a callout only when a
reader needs a distinct instruction to use an example correctly.

## Build and review

Install the documentation tools and run the checks from the repository root:

```sh
python -m pip install -e '.[dev,docs]'
./docs.sh
python scripts/check_docs_examples.py
./lint.sh
./test.sh
```

The example check requires human release 93 installed. `./docs.sh` builds
strictly and validates rendered internal links, anchors, README links and API
coverage. CI performs the docs build on PRs and publishes the site from main.
Use `mkdocs serve` to preview the pages; review the first-use and complete
reference pages at desktop and narrow widths, including long signatures and
tables. Do not use tests that only repeat documentation configuration values.

# Test design lessons

- Keep release and archive tests deterministic and offline. Exercise the real parser on a tiny local GTF, and use independent published filename examples plus a simulated platform listing to cover URL selection.

- Existing user caches are a compatibility contract. A cache-key fix must keep
  installed files and derived indexes usable in place, including offline and
  read-only/shared setups. Do not require reinstalling, renaming, migrating,
  rewriting, or rebuilding existing caches unless explicitly requested. Test
  preservation of paths, contents, modes, inodes, and timestamps before release.
  Apply new identity rules to new entries; document legacy identity limits
  instead of silently changing the assumptions of an existing installation.

- Distinguish missing registry coverage from information a caller must supply.
  Human Ensembl 75 already selects GRCh37; a dated selection with a known
  assembly name can resolve its accession from verified metadata. Never imply
  that all older assemblies require callers to know an accession. Treat
  explicitly supplied assembly/species/accession values as consistency
  constraints: reject known conflicts instead of silently overriding them.
  Distinguish proposed validation from checks the shipped API actually performs,
  and do not equate annotation dates with numbered releases or merge their
  caches without verifying dataset equivalence.

- Make assembly inference the normal interface: resolve numbered releases from
  species/release metadata and dated annotations from verified dataset records.
  Request an assembly only for genuine ambiguity, a nondefault selection, or
  an unregistered dataset. Accept broad names, exact patched names and versioned
  accessions as selectors; normalize their metadata separately. A broad name
  such as GRCh37 does not identify one patch/accession across every release.
  Explain which object attributes exist today, which are proposed, and where
  unknown legacy/custom metadata must stay None to preserve offline use.

- When the user explicitly says an unused compatibility API can be removed,
  retire its class/module/export rather than keeping a shim by default. Treat
  that as permission to remove its import and serialized-class compatibility,
  while retaining separately requested cache-file and active-API guarantees.
  Exercise argument combinations by installing and querying real tiny caches;
  checking constructor paths alone does not prove correct cache reuse.

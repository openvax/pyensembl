# Test design lessons

- Keep release and archive tests deterministic and offline. Exercise the real parser on a tiny local GTF, and use independent published filename examples plus a simulated platform listing to cover URL selection.

- Existing user caches are a compatibility contract. A cache-key fix must keep
  installed files and derived indexes usable in place, including offline and
  read-only/shared setups. Do not require reinstalling, renaming, migrating,
  rewriting, or rebuilding existing caches unless explicitly requested. Test
  preservation of paths, contents, modes, inodes, and timestamps before release.
  Apply new identity rules to new entries; document legacy identity limits
  instead of silently changing the assumptions of an existing installation.

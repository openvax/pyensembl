"""Custom fixtures must leave neither persistent caches nor source indexes."""

import json
import os
from pathlib import Path
import subprocess
import sys

from .data import data_path


def test_custom_fixtures_use_temporary_caches_and_clean_up(tmp_path):
    source_directory = Path(data_path(""))
    before = {path.name: path.stat().st_mtime_ns for path in source_directory.iterdir()}
    user_cache = tmp_path / "user-cache"
    result = subprocess.run(
        [sys.executable, "-c", """
import json
from pathlib import Path
from tests.data import custom_mouse_genome_grcm38_subset as mouse
from tests.test_tair10_complete import custom_tair10_genome_subset as tair
from tests.test_versions import mouse_genome, tair_genome

caches = []
for genome in (mouse, tair, mouse_genome, tair_genome):
    genome.index()
    assert genome.genes()
    cache = Path(genome.download_cache.cache_directory_path)
    assert Path(genome.db.local_db_path).parent == cache
    assert all(Path(p).parent == cache for p in
               genome.transcript_sequences.fasta_dictionary_pickle_paths)
    caches.append(str(cache))
    genome.close()
print(json.dumps(caches))
"""],
        env={**os.environ, "PYENSEMBL_CACHE_DIR": str(user_cache)},
        capture_output=True,
        text=True,
        check=True,
    )
    caches = json.loads(result.stdout)
    assert len(set(caches)) == 4
    assert not any(Path(path).exists() for path in caches)
    assert not user_cache.exists()
    after = {path.name: path.stat().st_mtime_ns for path in source_directory.iterdir()}
    assert after == before

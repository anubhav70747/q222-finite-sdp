"""Verify preserved package bytes and archives; not a mathematical proof audit."""
from pathlib import Path
import hashlib
import json
import zipfile

root = Path(__file__).resolve().parent
receipt = json.loads((root / 'IMPORT_VERIFICATION.json').read_text(encoding='utf-8'))
for row in receipt['files']:
    data = (root / row['destination']).read_bytes()
    assert len(data) == row['bytes'], row['destination']
    assert hashlib.sha256(data).hexdigest() == row['sha256'], row['destination']
members = 0
for row in receipt['archives']:
    with zipfile.ZipFile(root / 'archives' / row['file']) as archive:
        assert archive.testzip() is None
        assert len(archive.namelist()) == row['members']
        for name in archive.namelist():
            destination = 'package-readme.txt' if name == 'README.txt' else name
            assert archive.read(name) == (root / destination).read_bytes(), name
            members += 1
print(json.dumps({'ok': True, 'files': len(receipt['files']),
                  'archive_members': members,
                  'scope': 'preservation only; no new mathematical certification'}))

#!/usr/bin/env python3
"""Verify every new claim and the reused validity identity, exactly.

Python 3 standard library only. No optimization, tolerance, CAS, or network.
The included certificate and matrix data are not modified.
"""
from pathlib import Path
import hashlib
import json
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
DATA = ROOT/'data'
sys.dont_write_bytecode = True

def run(script, *args):
    result = subprocess.run([sys.executable, '-B', str(script), *map(str, args)],
                            cwd=ROOT, text=True, capture_output=True)
    if result.returncode:
        print(result.stdout)
        print(result.stderr, file=sys.stderr)
        raise SystemExit(result.returncode)
    print('PASS:', script.name, flush=True)
    return result.stdout

run(DATA/'verify_exact_separator.py', DATA/'dyadic_order33_separator.json')
run(DATA/'verify_ramified.py')
run(DATA/'verify_ramified_pencil.py')
run(DATA/'verify_ramified_pencil_independent.py')
run(DATA/'build_function_and_equation_data.py', '--check')

# Rebuild both 148-entry auxiliary encodings in a private temporary copy,
# then compare the exact JSON against the supplied evaluated values.
with tempfile.TemporaryDirectory(prefix='q222-halfline-') as tmp:
    work = Path(tmp)
    for name in ['extract_ramified_auxiliaries.py', 'ramified_order49_pencil.json',
                 'ramified_exact.json', 'ramified_cbar_exact.json']:
        (work/name).write_bytes((DATA/name).read_bytes())
    for source, output in [('ramified_exact.json', 'ramified_order49_source_auxiliaries.json'),
                           ('ramified_cbar_exact.json', 'ramified_order49_cbar_auxiliaries.json')]:
        run(work/'extract_ramified_auxiliaries.py', work/source)
        assert json.loads((work/output).read_text()) == json.loads((DATA/output).read_text())
    # The published small source certificate is retained as an established
    # input. Its checker can also verify the determinant directly here.
    source = ROOT/'reused/ESTABLISHED_VALID_SOURCE.json'
    run(ROOT/'reused/check_source_validity.py', '--source', source,
        '--public-source', source, '--output-dir', work)
    generated = json.loads((work/'pro_completion_selfcontained_data.json').read_text())
    retained = json.loads((ROOT/'reused/pro_completion_selfcontained_data.json').read_text())
    assert generated == retained

manifest_path = ROOT/'SHA256SUMS.json'
if manifest_path.exists():
    for name, digest in json.loads(manifest_path.read_text()).items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == digest, name
    print('PASS: supplied file hashes', flush=True)

print('RESULT: C33 is strictly contained in the ramified order-49 cone,')
print('which is contained in the full nonnegative canonical valid-source cone.')
print('Global adequacy of the ramified cone remains unproved.')

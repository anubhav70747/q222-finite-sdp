#!/usr/bin/env python3
"""Reproduce the new exact checks; optionally rebuild only the new matrices."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parent
parser=argparse.ArgumentParser()
parser.add_argument('--rebuild',action='store_true')
args=parser.parse_args()

def run(relative,*arguments):
    path=ROOT/relative
    print('\nRUN '+str(relative),flush=True)
    subprocess.run([sys.executable,str(path),*arguments],cwd=path.parent,check=True)

p=json.loads((ROOT/'PROVENANCE.json').read_text())
old=ROOT/'reused/C1_UNIFORM_LOCAL_ORDER669.coo.json'
assert hashlib.sha256(old.read_bytes()).hexdigest()==p['reused_order669_sha256']
print('PASS: unchanged order-669 matrix SHA256; no old compilation run.',flush=True)
run('weighted/verify_global_multiplier.py')
run('weighted/check_determinant_data.py')
if args.rebuild:
    run('weighted/compile_weighted_family.py')
    run('weighted/compile_weighted_pencil.py')
run('weighted/check_evaluated_matrices.py','--full-source-identity')
if args.rebuild:
    run('angular-order49/compile_angular_pencil.py')
run('angular-order49/check_angular_pencil.py')
run('bridges/bounded_R/check_bridge.py')
run('bridges/compact_Bob/check_compact_bob.py')
run('bridges/bounded_unit_L/check_bounded_unit_L.py')
run('bridges/bounded_unit_L/check_independent_annulus.py')
print('\nAll new exact checks passed. The full affine Q222 lift remains open.',flush=True)

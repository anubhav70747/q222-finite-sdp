#!/usr/bin/env python3
"""Run the exact coefficient/Schur/native replay. No search or numerical SDP."""
import subprocess,sys,pathlib
root=pathlib.Path(__file__).resolve().parent
steps=[
 'verify_z18_exact.py',
 'verify_schur_identities.py',
 'original_baseline/verify_c1_native_stationarity.py',
 'original_baseline/verify_c1_map_exact.py',
 'original_baseline/verify_parent_scale.py',
]
for step in steps:
 p=subprocess.run([sys.executable,str(root/step)],cwd=root,capture_output=True,text=True)
 print('===',step,'exit',p.returncode,'===',flush=True)
 for line in p.stdout.splitlines():
  if ('PASS' in line or 'EXACT_' in line or 'ALL ' in line or 'rank ' in line or 'No universal' in line):print(line,flush=True)
 if p.returncode:
  print(p.stdout,p.stderr,flush=True)
  raise SystemExit(p.returncode)
print('EXACT REPLAY SUCCESS. Universal oriented-source PSD is NOT established.',flush=True)

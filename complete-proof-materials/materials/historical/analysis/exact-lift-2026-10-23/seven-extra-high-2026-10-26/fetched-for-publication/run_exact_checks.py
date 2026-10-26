#!/usr/bin/env python3
"""Deterministic exact-rational/polynomial replay, no sampling or SDP solver."""
import subprocess,sys,pathlib
here=pathlib.Path(__file__).resolve().parent
names=['Q222_R6_R7_source_exact_audit.py',
       'Q222_R6R7_full_parent_readback_audit.py',
       'Q222_R6R7_R34_identity_audit.py',
       'Q222_R34_short_balance_audit.py',
       'Q222_R34_full_parent_readback_audit.py']
for name in names:
    print('>>>',name,flush=True)
    subprocess.run([sys.executable,str(here/name)],check=True,timeout=90)
print('ALL FIVE EXACT AUDITS PASSED',flush=True)

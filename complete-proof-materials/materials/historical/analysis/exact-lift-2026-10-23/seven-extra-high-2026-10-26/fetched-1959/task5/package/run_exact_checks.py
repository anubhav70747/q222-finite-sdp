#!/usr/bin/env python3
"""Run only the three deterministic symbolic proof identities, record a reproducible log."""
from pathlib import Path
import subprocess,sys
root=Path(__file__).resolve().parent
lines=[]
for script in ('verify_full_parent_bob.py','verify_boundary_schur.py','verify_anchor_interpolation.py'):
    p=subprocess.run([sys.executable,str(root/script)],cwd=str(root),capture_output=True,text=True,timeout=90)
    lines.append('=== '+script+' ===\n'+p.stdout+p.stderr+'exit='+str(p.returncode)+'\n')
    if p.returncode:
        (root/'verification.log').write_text('\n'.join(lines),encoding='utf-8')
        raise RuntimeError(script+' failed')
(root/'verification.log').write_text('\n'.join(lines),encoding='utf-8')
print('PASS three exact full-state symbolic modules; log written.')

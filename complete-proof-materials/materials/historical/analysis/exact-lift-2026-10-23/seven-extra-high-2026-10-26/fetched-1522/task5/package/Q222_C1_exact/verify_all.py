#!/usr/bin/env python3
"""Focused exact identity-only replay. No solver, samples, or searches."""
import subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent
for script in ['verify_parent_scale.py','verify_c1_native_stationarity.py','verify_c1_map_exact.py','verify_forced_positivity.py']:
 print('\n===== '+script+' =====',flush=True)
 subprocess.run([sys.executable,str(ROOT/script)],cwd=str(ROOT),check=True)
print('\nALL EXACT CHECKERS PASS. THIS DOES NOT ESTABLISH ANY FULL 7x7 OR 5x5 PSD COMPLETION.',flush=True)

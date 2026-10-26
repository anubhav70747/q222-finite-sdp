"""Run the three exact Q222 C1 BN/UN checkers. Zero exit status means each claim checked algebraically."""
from pathlib import Path
import runpy
here=Path(__file__).resolve().parent
for name in ['verify_bn_un_independent.py','verify_bn_un_parent_full.py','verify_bn_un_structural_coeffs.py']:
 print('\n===',name,'===', flush=True)
 runpy.run_path(str(here/name),run_name='__main__')
print('\nALL THREE BN/UN EXACT CHECKERS PASSED.')

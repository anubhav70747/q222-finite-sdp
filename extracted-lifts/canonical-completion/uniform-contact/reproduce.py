"""Reproduce the complete fixed-cap certificate and its rational matrices."""
from pathlib import Path
import subprocess,sys
HERE=Path(__file__).resolve().parent
for script in [
    'check_source_parent.py',
    'check_normal_and_ru_budgets.py',
    'extraction/build_support.py',
    'extraction/check_factored_assembly.py',
    'extraction/compile_affine.py',
    'extraction/check_fixed_grid.py',
    'check_annulus_obstruction.py',
]:
    print('\n'+script,flush=True)
    subprocess.run([sys.executable,'-B',str(HERE/script)],cwd=HERE,check=True)
print('\nAll exact checks completed. Scope: the full-R1 fixed testing cap.')

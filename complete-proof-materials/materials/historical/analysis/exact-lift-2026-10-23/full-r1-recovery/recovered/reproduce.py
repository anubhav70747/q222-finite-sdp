#!/usr/bin/env python3
"""Reproduce only this new all-R1 cap construction (Python + SymPy)."""
from pathlib import Path
import subprocess
import sys

root=Path(__file__).resolve().parent
steps=[('build_cap_family.py',[]),('compile_cap_affine.py',[]),
       ('r1_sos_roots.py',[]),('check_complete_budgets.py',[]),
       ('check_uniform_normal_core.py',[]),
       ('review_full_parent_recombination.py',[]),
       ('check_fixed_cap_support.py',[]),('check_cap_closure_grid.py',[]),
       ('check_cap_closure_grid.py',['--matrices'])]
(root/'data').mkdir(exist_ok=True)
log=[]
for name,args in steps:
    run=subprocess.run([sys.executable,str(root/'scripts'/name),*args],
                       capture_output=True,text=True)
    print(run.stdout,end='')
    if run.returncode:
        print(run.stderr,file=sys.stderr)
        raise SystemExit(run.returncode)
    log.append(name+' '+' '.join(args)+'\n'+run.stdout)
(root/'data'/'verification.txt').write_text('\n'.join(log))
if '--pdf' in sys.argv:
    for _ in range(2):
        subprocess.run(['pdflatex','-interaction=nonstopmode','-halt-on-error',
                        'all_r1_cap.tex'],cwd=root/'proof',check=True,
                       stdout=subprocess.DEVNULL)
print('Complete: fixed local cap family and exact source-affine records.')
print('Global canonical coverage has not been established.')

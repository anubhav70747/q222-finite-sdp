"""Run the fresh diagnostic programs; no archived program is imported."""
from pathlib import Path
import shutil, subprocess, sys
import sympy

if sympy.__version__!='1.14.0':
    raise SystemExit('Use SymPy 1.14.0: python -m pip install sympy==1.14.0')
root=Path(__file__).resolve().parent
dest=root/'recomputed';dest.mkdir(exist_ok=True)
names=['root_identity_audit','root_boundary_identity_audit','root_native_entry_audit',
       'positive_scale_independent','strong_c2_independent','unbounded_b_identity_audit',
       'local_identity_audit','support_truth_tables']
for name in names:
    p=dest/(name+'.py');shutil.copy2(root/'programs'/p.name,p)
    cmd=[sys.executable,str(p)]
    if name=='local_identity_audit':cmd+=['--all','--log-json',str(dest/(name+'.checks.json'))]
    with (dest/(name+'.log')).open('w') as log:
        result=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,cwd=dest)
    if result.returncode:
        raise SystemExit('Failed: '+name+'; inspect '+str(dest/(name+'.log')))
    print('PASS '+name,flush=True)
print('Completed. Original receipts are unchanged; fresh outputs are in checks/recomputed/.')

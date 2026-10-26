"""Sequential exact rational/symbolic identity replay; no solvers/sweeps."""
from pathlib import Path
import subprocess, sys, time
names=[
 'Q222_IC_LII_small_identity_check.py',
 'Q222_IC_normal_symbolic_check.py',
 'Q222_IC_graph_support_check.py',
 'Q222_IH_endpoint_source_check.py',
 'Q222_WIE4_full_graph_identity.py',
 'Q222_LII_full_parent_readback.py',
 'Q222_LII_coeff_extraction_probe.py',
 'Q222_IC_LII_residue_check.py',
]
root=Path(__file__).resolve().parent
for f in names:
    start=time.monotonic()
    proc=subprocess.run([sys.executable,str(root/f)],capture_output=True,text=True,timeout=40)
    if proc.returncode:
        print(proc.stdout); print(proc.stderr,file=sys.stderr)
        raise SystemExit(f'FAILED: {f}')
    print(f'PASS {f} ({time.monotonic()-start:.1f}s)',flush=True)
print('ALL EXACT CHECKS PASSED')

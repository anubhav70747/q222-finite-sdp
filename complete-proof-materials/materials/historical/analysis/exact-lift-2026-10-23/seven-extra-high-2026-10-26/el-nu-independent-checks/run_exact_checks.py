import subprocess,sys,pathlib
p=pathlib.Path(__file__).resolve().parent
for f in ['verify_sd_source.py','verify_sphere.py','verify_parent.py','verify_el_parent_verified.py','check_el20.py','verify_budgets.py','verify_ce.py','verify_eb_fu.py']:
 print('##',f,flush=True)
 r=subprocess.run([sys.executable,str(p/f)],capture_output=True,text=True,timeout=65)
 print(r.stdout,flush=True)
 if r.returncode!=0 or any(x in r.stdout for x in ['FAIL','False']):
  print('ERROR in',f, r.stderr)
  sys.exit(1)
print('ALL EXACT CHECKS PASSED')

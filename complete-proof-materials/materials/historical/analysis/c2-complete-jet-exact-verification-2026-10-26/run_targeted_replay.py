from pathlib import Path
import json,zipfile,hashlib,shutil,subprocess,sys,os,time,ast
ROOT=Path(r'C:/Users/anubh/OneDrive/Documents/ChatGPT/ManuscriptWriting/q222/analytical-paper-2026-10-22')
SRC=ROOT/'analysis/exact-lift-2026-10-23/seven-extra-high-2026-10-26/fetched-1939/task3'
OUT=Path(__file__).resolve().parent
import sympy
r={'python':sys.version,'sympy':sympy.__version__,'environment':{}}
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','PYTHONDONTWRITEBYTECODE','PYTHONUTF8'):
    os.environ[key]='1';r['environment'][key]='1'
arc=SRC/'original-package.zip';raw=arc.read_bytes();r['archive_bytes']=len(raw);r['archive_sha256']=hashlib.sha256(raw).hexdigest()
with zipfile.ZipFile(arc) as z:
    r['crc_bad_member']=z.testzip();files=[x for x in z.infolist() if not x.is_dir()];r['archive_file_count']=len(files)
    r['extracted_byte_checks']=[{'file':x.filename,'bytes':x.file_size,'identical':(SRC/'package'/x.filename).read_bytes()==z.read(x)} for x in files]
m=json.loads((SRC/'package/MANIFEST.json').read_text())
r['manifest_checks']=[{'file':name,'sha256_match':hashlib.sha256((SRC/'package'/name).read_bytes()).hexdigest()==spec['sha256'],'bytes_match':(SRC/'package'/name).stat().st_size==spec['bytes']} for name,spec in m['files'].items()]
assert r['crc_bad_member'] is None and all(x['identical'] for x in r['extracted_byte_checks']) and all(x['sha256_match'] and x['bytes_match'] for x in r['manifest_checks'])
print('PASS archive CRC, all 4 extracted bytes and all 3 manifest SHA256/byte checks',flush=True)
private=OUT/'private-copy'
if private.exists():raise RuntimeError('Private copy exists')
shutil.copytree(SRC/'package',private)
r['script_assert_statement_count']=sum(isinstance(n,ast.Assert) for n in ast.walk(ast.parse((private/'verify_exact.py').read_text())))
start=time.perf_counter()
with (OUT/'verify_exact.log').open('w',encoding='utf-8') as f:
    p=subprocess.run([sys.executable,'-B','-u','verify_exact.py'],cwd=private,env=os.environ.copy(),stdout=f,stderr=subprocess.STDOUT,timeout=600)
r['replay']={'exit_code':p.returncode,'elapsed_seconds':round(time.perf_counter()-start,6)}
print('REPLAY',r['replay'],flush=True)
if p.returncode:print((OUT/'verify_exact.log').read_text(),flush=True)
assert p.returncode==0
original=(SRC/'package/EXACT_RESULTS.json').read_bytes();fresh=(private/'EXACT_RESULTS.json').read_bytes()
r['generated_data_comparison']={'byte_identical':original==fresh,'newline_normalized_identical':original.replace(b'\r\n',b'\n')==fresh.replace(b'\r\n',b'\n'),'json_semantically_identical':json.loads(original)==json.loads(fresh)}
r['reported_pass_lines']=sum(line.startswith('PASS:') for line in (OUT/'verify_exact.log').read_text().splitlines())
assert r['generated_data_comparison']['json_semantically_identical']
(OUT/'integrity_and_replay_results.json').write_text(json.dumps(r,indent=2),encoding='utf-8')
print('DATA',r['generated_data_comparison'],'PASS lines',r['reported_pass_lines'],'source assertions',r['script_assert_statement_count'],flush=True)

from pathlib import Path
import json,zipfile,hashlib,shutil,subprocess,sys,os,time,ast
ROOT=Path(r'C:/Users/anubh/OneDrive/Documents/ChatGPT/ManuscriptWriting/q222/analytical-paper-2026-10-22')
SRC=ROOT/'analysis/exact-lift-2026-10-23/seven-extra-high-2026-10-26/fetched-20261027-0705'
OUT=Path(__file__).resolve().parent
import sympy
r={'python':sys.version,'sympy':sympy.__version__,'environment':{},'integrity':{},'runs':[]}
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','PYTHONDONTWRITEBYTECODE','PYTHONUTF8'):
    os.environ[key]='1';r['environment'][key]='1'
names={'task5':'Q222_C1_K364_CONTACT_JET_EXACT_LOCAL_POSITIVITY_2026-10-27.zip','task7':'K364_MOMENT_DUAL_ANALYTIC_AUDIT.zip'}
def save():(OUT/'results.json').write_text(json.dumps(r,indent=2),encoding='utf-8')
for task in ('task5','task7'):
    arc=SRC/task/'original-package.zip';raw=arc.read_bytes();pkg=SRC/task/'package'
    zrec={'archive_bytes':len(raw),'archive_sha256':hashlib.sha256(raw).hexdigest(),'download_original':str(Path(r'C:/Users/anubh/Downloads')/names[task])}
    zrec['download_original_bytes_match']=Path(zrec['download_original']).read_bytes()==raw
    with zipfile.ZipFile(arc) as z:
        zrec['crc_bad_member']=z.testzip();members=[x for x in z.infolist() if not x.is_dir()]
        zrec['archive_file_count']=len(members)
        zrec['member_checks']=[{'file':info.filename,'bytes':info.file_size,'extracted_byte_identical':(pkg/info.filename).read_bytes()==z.read(info)} for info in members]
        zrec['pycache_members']=[info.filename for info in members if '__pycache__' in info.filename or info.filename.endswith('.pyc')]
    base=pkg/'Q222_C1_K364_CONTACT_JET_2026_10_27' if task=='task5' else pkg
    manifest=json.loads((base/('MANIFEST.json' if task=='task5' else 'SHA256_MANIFEST.json')).read_text())
    checks=[]
    for name,spec in (manifest['files'] if task=='task5' else manifest).items():
        data=(base/name).read_bytes()
        checks.append({'file':name,'sha256_match':hashlib.sha256(data).hexdigest()==(spec['sha256'] if task=='task5' else spec),'bytes_match':len(data)==spec['bytes'] if task=='task5' else None})
    zrec['manifest_checks']=checks
    assert zrec['download_original_bytes_match'] and zrec['crc_bad_member'] is None and all(x['extracted_byte_identical'] for x in zrec['member_checks']) and all(x['sha256_match'] and x['bytes_match'] is not False for x in checks)
    r['integrity'][task]=zrec;save()
    dest=OUT/'private-copy'/task
    if dest.exists():raise RuntimeError('Private copy exists')
    shutil.copytree(pkg,dest,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
    print(task,'PASS integrity',len(members),'members',len(checks),'manifest entries',zrec['archive_bytes'],'bytes',flush=True)
commands=[('task5/Q222_C1_K364_CONTACT_JET_2026_10_27','verify_contact_geometry.py'),('task5/Q222_C1_K364_CONTACT_JET_2026_10_27','verify_two_schur.py'),('task5/Q222_C1_K364_CONTACT_JET_2026_10_27','verify_contact_jet_exact.py'),('task7','verify_exact_dual.py')]
for folder,script in commands:
    cwd=OUT/'private-copy'/folder;log=OUT/(Path(script).stem+'.log');start=time.perf_counter()
    print('START',script,flush=True)
    with log.open('w',encoding='utf-8') as f:
        p=subprocess.run([sys.executable,'-B','-u',script],cwd=cwd,env=os.environ.copy(),stdout=f,stderr=subprocess.STDOUT,timeout=300)
    text=log.read_text();rec={'script':script,'exit_code':p.returncode,'elapsed_seconds':round(time.perf_counter()-start,6),'log':log.name,'pass_lines':sum(line.startswith('PASS') for line in text.splitlines()),'assert_statements':sum(isinstance(n,ast.Assert) for n in ast.walk(ast.parse((cwd/script).read_text())))}
    r['runs'].append(rec);save();print('END',rec,flush=True)
    if p.returncode:print(text,flush=True)
    assert p.returncode==0
old=(SRC/'task7/package/EXACT_EVIDENCE.json').read_bytes();new=(OUT/'private-copy/task7/EXACT_EVIDENCE.json').read_bytes()
r['task7_machine_comparison']={'byte_identical':old==new,'newline_normalized_identical':old.replace(b'\r\n',b'\n')==new.replace(b'\r\n',b'\n'),'json_semantically_identical':json.loads(old)==json.loads(new)}
assert r['task7_machine_comparison']['json_semantically_identical']
original_log=(SRC/'task5/package/Q222_C1_K364_CONTACT_JET_2026_10_27/verification.log').read_text()
fresh_log=''.join((OUT/(Path(script).stem+'.log')).read_text() for folder,script in commands if folder.startswith('task5'))
r['task5_combined_log_exact_text_match']=original_log==fresh_log
save();print('MACHINE COMPARISONS',r['task7_machine_comparison'],'task5 combined log match',r['task5_combined_log_exact_text_match'],flush=True)

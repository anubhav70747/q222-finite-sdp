from pathlib import Path
import hashlib, json, zipfile, shutil, sys, os, subprocess, time
ROOT=Path(r'C:/Users/anubh/OneDrive/Documents/ChatGPT/ManuscriptWriting/q222/analytical-paper-2026-10-22')
SRC=ROOT/'analysis/exact-lift-2026-10-23/seven-extra-high-2026-10-26/fetched-1739'
OUT=ROOT/'analysis/latest-transfer-exact-verification-2026-10-26'
report={'python':sys.version,'environment':{},'integrity':{},'replays':[]}
import sympy
report['sympy']=sympy.__version__
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS'):
    os.environ[key]='1';report['environment'][key]='1'
for task in ('task3','task5'):
    arc=SRC/task/'original-package.zip';pkg=SRC/task/'package'
    result={'archive':str(arc),'archive_bytes':arc.stat().st_size,'archive_sha256':hashlib.sha256(arc.read_bytes()).hexdigest()}
    with zipfile.ZipFile(arc) as z:
        result['crc_bad_member']=z.testzip()
        members=[x for x in z.infolist() if not x.is_dir()]
        mismatches=[]
        for info in members:
            p=pkg/info.filename
            if not p.is_file():mismatches.append({'member':info.filename,'problem':'missing extracted file'})
            elif p.read_bytes()!=z.read(info):mismatches.append({'member':info.filename,'problem':'extracted bytes mismatch'})
        result['archive_file_count']=len(members);result['extracted_byte_mismatches']=mismatches
    base=pkg/'C2_Strict_Coupled_Transfer_20261026' if task=='task3' else pkg
    manifest=base/('MANIFEST.json' if task=='task3' else 'HASHES.json')
    m=json.loads(manifest.read_text(encoding='utf-8'))
    entries=[{'path':p,**item} for p,item in m['files'].items()] if task=='task3' else m['files']
    checks=[]
    for entry in entries:
        p=base/entry['path'];data=p.read_bytes() if p.is_file() else None
        checks.append({'path':entry['path'],'exists':data is not None,'bytes_match':data is not None and len(data)==entry['bytes'],'sha256_match':data is not None and hashlib.sha256(data).hexdigest()==entry['sha256']})
    result['manifest_checks']=checks
    result['all_integrity_ok']=result['crc_bad_member'] is None and not mismatches and all(x['exists'] and x['bytes_match'] and x['sha256_match'] for x in checks)
    report['integrity'][task]=result
    print(task,'archive',result['archive_bytes'],'bytes;',len(members),'files;',len(checks),'manifest entries; integrity',result['all_integrity_ok'],flush=True)
    assert result['all_integrity_ok'],result
    private=OUT/'private-copy'/task
    if private.exists():raise RuntimeError('Private copy already exists: '+str(private))
    shutil.copytree(pkg,private)
(OUT/'integrity.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
commands=[('C2','task3/C2_Strict_Coupled_Transfer_20261026','verify_c2_transfer.py'),('C1_Z18','task5','verify_z18_exact.py'),('C1_Schur','task5','verify_schur_identities.py')]
for label,sub,script in commands:
    cwd=OUT/'private-copy'/sub
    log=OUT/(label+'.log');start=time.monotonic()
    print('START',label,flush=True)
    with log.open('w',encoding='utf-8') as f:
        proc=subprocess.run([sys.executable,'-u',script],cwd=cwd,env=os.environ.copy(),stdout=f,stderr=subprocess.STDOUT,timeout=600)
    result={'label':label,'script':script,'cwd':str(cwd),'log':str(log),'exit_code':proc.returncode,'elapsed_seconds':round(time.monotonic()-start,3)}
    report['replays'].append(result)
    (OUT/'results.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print('END',label,result,flush=True)
print('TARGETED REPLAY COMPLETE',flush=True)

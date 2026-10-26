from pathlib import Path
import hashlib, json, zipfile, shutil, sys, os, subprocess, time, gzip
ROOT=Path(r'C:/Users/anubh/OneDrive/Documents/ChatGPT/ManuscriptWriting/q222/analytical-paper-2026-10-22')
SRC=ROOT/'analysis/exact-lift-2026-10-23/seven-extra-high-2026-10-26/fetched-1819'
OUT=ROOT/'analysis/global-pencil-exact-verification-2026-10-26'
import sympy
report={'python':sys.version,'sympy':sympy.__version__,'environment':{},'integrity':{},'replays':[],'comparisons':[]}
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','PYTHONDONTWRITEBYTECODE','PYTHONUTF8'):
    os.environ[key]='1';report['environment'][key]='1'
receipts={f"task{x['role']}":x for x in json.loads((SRC/'RECOVERY_RECEIPT.json').read_text())}
def save():
    (OUT/'results.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
for task in ('task7','task5'):
    arc=SRC/task/'original-package.zip';pkg=SRC/task/'package';raw=arc.read_bytes();sha=hashlib.sha256(raw).hexdigest()
    result={'archive':str(arc),'archive_bytes':len(raw),'archive_sha256':sha,'receipt_sha256_match':sha==receipts[task]['sha256']}
    original=Path(receipts[task]['source'])
    result['download_original']=str(original)
    result['download_original_bytes_identical']=original.read_bytes()==raw
    with zipfile.ZipFile(arc) as z:
        result['crc_bad_member']=z.testzip();members=[x for x in z.infolist() if not x.is_dir()]
        mismatches=[]
        for info in members:
            p=pkg/info.filename
            if not p.is_file():mismatches.append({'member':info.filename,'problem':'missing extracted file'})
            elif p.read_bytes()!=z.read(info):mismatches.append({'member':info.filename,'problem':'extracted bytes mismatch'})
        result['archive_file_count']=len(members);result['extracted_byte_mismatches']=mismatches
        result['pycache_archived_members']=[x.filename for x in members if '__pycache__' in x.filename or x.filename.endswith('.pyc')]
    base=pkg if task=='task7' else pkg/'Q222_C1_ORIENTED_K364_ADMISSION_2026_10_26'
    if task=='task7':
        entries=[{'path':line.split(None,1)[1].strip(),'sha256':line.split(None,1)[0]} for line in (base/'MANIFEST.sha256').read_text().splitlines() if line.strip()]
    else:entries=json.loads((base/'MANIFEST.json').read_text())['entries']
    checks=[]
    for entry in entries:
        p=base/entry['path'];data=p.read_bytes() if p.is_file() else None
        checks.append({'path':entry['path'],'exists':data is not None,'bytes':len(data) if data is not None else None,'bytes_match':len(data)==entry['bytes'] if data is not None and 'bytes' in entry else None,'sha256_match':data is not None and hashlib.sha256(data).hexdigest()==entry['sha256'],'executable_evidence_excluded':'__pycache__' in entry['path'] or entry['path'].endswith('.pyc')})
    result['manifest_checks']=checks
    result['all_integrity_ok']=result['receipt_sha256_match'] and result['download_original_bytes_identical'] and result['crc_bad_member'] is None and not mismatches and all(x['exists'] and x['bytes_match'] is not False and x['sha256_match'] for x in checks)
    report['integrity'][task]=result;save()
    print(task,'archive',result['archive_bytes'],'bytes;',len(members),'files;',len(checks),'manifest entries; integrity',result['all_integrity_ok'],flush=True)
    assert result['all_integrity_ok'],result
    private=OUT/'private-copy'/task
    if private.exists():raise RuntimeError('Private copy already exists: '+str(private))
    shutil.copytree(pkg,private,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
def run(label,sub,script,args=()):
    cwd=OUT/'private-copy'/sub;log=OUT/(label+'.log');start=time.monotonic()
    print('START',label,flush=True)
    with log.open('w',encoding='utf-8') as f:
        proc=subprocess.run([sys.executable,'-B','-u',script,*args],cwd=cwd,env=os.environ.copy(),stdout=f,stderr=subprocess.STDOUT,timeout=600)
    result={'label':label,'script':script,'args':list(args),'cwd':str(cwd),'log':str(log),'exit_code':proc.returncode,'elapsed_seconds':round(time.monotonic()-start,3)}
    report['replays'].append(result);save();print('END',label,result,flush=True)
    if proc.returncode:print(log.read_text(),flush=True);raise RuntimeError(label+' failed')
run('C1_source_charge','task5/Q222_C1_ORIENTED_K364_ADMISSION_2026_10_26','verify_source_charge_margin.py')
run('global_exact_checker','task7','verify_exact_pencil.py')
# Regenerate once, in the same disposable copy, after preserving fresh checker outputs.
p7=OUT/'private-copy/task7'
for name in ('EXACT_VERIFICATION_REPORT.json','K364_STRICT_RATIONAL_SOURCE_AND_GRAM.json'):
    shutil.copy2(p7/name,OUT/('fresh_'+name))
run('global_generate','task7','generate_exact_pencil.py')
run('global_assemble','task7','assemble_global_pencil.py',('--check-structural',))
run('global_split','task7','split_materialized_tables.py')
for orig in sorted((SRC/'task7/package').glob('K364*.json*')):
    now=p7/orig.name;oldbytes=orig.read_bytes();newbytes=now.read_bytes()
    entry={'file':orig.name,'original_bytes':len(oldbytes),'generated_bytes':len(newbytes),'byte_identical':oldbytes==newbytes,'original_sha256':hashlib.sha256(oldbytes).hexdigest(),'generated_sha256':hashlib.sha256(newbytes).hexdigest()}
    oldjson=gzip.decompress(oldbytes) if orig.suffix=='.gz' else oldbytes
    newjson=gzip.decompress(newbytes) if now.suffix=='.gz' else newbytes
    entry['decompressed_bytes_identical']=oldjson==newjson
    entry['json_semantically_identical']=json.loads(oldjson)==json.loads(newjson)
    entry['newline_normalized_identical']=oldjson.replace(b'\r\n',b'\n')==newjson.replace(b'\r\n',b'\n')
    report['comparisons'].append(entry)
    print('COMPARE',orig.name,'bytes',entry['byte_identical'],'JSON',entry['json_semantically_identical'],flush=True)
    assert entry['json_semantically_identical'],entry
save()
print('TARGETED CHECKS AND ONE REGENERATION COMPLETE',flush=True)

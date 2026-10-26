#!/usr/bin/env python3
"""Verify the delivered files; optionally replay exact suites in isolation."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import zipfile

ROOT=Path(__file__).resolve().parent

def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()

def integrity(root,manifest):
    mismatches=[]
    for row in manifest['files']:
        p=root/row['path']
        if not p.is_file():mismatches.append({'path':row['path'],'error':'missing'})
        elif p.stat().st_size!=row['bytes'] or digest(p)!=row['sha256']:
            mismatches.append({'path':row['path'],'error':'content mismatch'})
    return mismatches

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--replay',choices=['quick'])
    parser.add_argument('--report',type=Path,help='Write this run result to a separate report path.')
    args=parser.parse_args()
    if not __debug__:
        raise SystemExit('Do not run verification with Python -O; exact assertions are required.')
    manifest=json.loads((ROOT/'MANIFEST.json').read_text())
    bad=integrity(ROOT,manifest)
    result={'manifest_status':'PASS' if not bad else 'FAIL','manifest_bound_files':len(manifest['files']),'mismatches':bad,'global_E1_E2_E3':'OPEN','historical_549_page_existence_proof_reaudited':False}
    if bad:
        print(json.dumps(result,indent=2,sort_keys=True));return 1
    inherited=ROOT/'recovered/Q222_Extraction_Continuation_Package_2026-10-22.zip'
    with zipfile.ZipFile(inherited) as z:
        err=z.testzip()
        if err:raise AssertionError('Retained archive CRC failure: '+err)
    if digest(inherited)!='9dcf9554b1a2578bb6b9fd909e7d194d82636fe61a0035da0694eb2e44f441a3':
        raise AssertionError('Retained archive source hash mismatch')
    result['retained_archive_integrity']='PASS'
    result['retained_archive_sha256']=digest(inherited)
    if args.replay:
        plan=json.loads((ROOT/'REPLAY_PLAN.json').read_text())
        rows=[]
        with tempfile.TemporaryDirectory(prefix='q222_regenerated_replay_') as tmp:
            target=Path(tmp)/'package'
            # Copy only manifest-bound deliverables, rather than scratch files.
            for row in manifest['files']:
                src=ROOT/row['path'];dst=target/row['path']
                dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dst)
            (target/'MANIFEST.json').write_bytes((ROOT/'MANIFEST.json').read_bytes())
            env=dict(os.environ);env['PYTHONDONTWRITEBYTECODE']='1'
            for suite in plan['suites']:
                proc=subprocess.run([sys.executable,str(target/suite['script'])],cwd=target,env=env,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=300)
                row={'name':suite['name'],'script':suite['script'],'exit_code':proc.returncode,'status':'PASS' if proc.returncode==0 else 'FAIL','stdout':proc.stdout.strip(),'stderr':proc.stderr.strip()}
                rows.append(row)
                if proc.returncode!=0:
                    result['replay']={'status':'FAIL','suites':rows}
                    if args.report:
                        args.report.parent.mkdir(parents=True,exist_ok=True);args.report.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
                    print(json.dumps(result,indent=2,sort_keys=True));return 1
            changed=integrity(target,manifest)
            result['replay']={'status':'PASS' if not changed else 'FAIL','suites':rows,'suite_count':len(rows),'manifest_bound_output_changes':changed,'isolation':'All suites run in a temporary copied package; delivered files are not rewritten.'}
        if changed:
            if args.report:
                args.report.parent.mkdir(parents=True,exist_ok=True);args.report.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
            print(json.dumps(result,indent=2,sort_keys=True));return 1
    result['status']='PASS'
    result['scope']='Integrity and the named finite exact checks; no global extraction completion, external peer review, or proof-assistant verification.'
    if args.report:
        args.report.parent.mkdir(parents=True,exist_ok=True);args.report.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps(result,indent=2,sort_keys=True))
    return 0

if __name__=='__main__':raise SystemExit(main())

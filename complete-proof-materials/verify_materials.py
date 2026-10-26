"""Verify published bytes, source identities and checklist links; no theorem checker."""
from pathlib import Path
import ast
import hashlib
import json
import re
import sys
import zipfile
from urllib.parse import unquote

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
manifest = json.loads((HERE/'MANIFEST.json').read_text(encoding='utf-8'))
errors = []
verified = 0
for entry in manifest['files']:
    p = REPO / entry['path']
    if not p.is_file():
        errors.append('missing: ' + entry['path'])
        continue
    data=p.read_bytes()
    if entry.get('canonical_lf'):
        data=data.replace(b'\r\n',b'\n')
    if len(data)!=entry['bytes'] or hashlib.sha256(data).hexdigest()!=entry['sha256']:
        errors.append('bytes changed: '+entry['path'])
    else:
        verified+=1

checklist=json.loads((HERE/'CHECKLIST.json').read_text(encoding='utf-8'))
for c in checklist['components']:
    for kind,field in [('source','source_sha256'),('presentation','presentation_sha256')]:
        p=HERE/'materials/historical/recovery/full_r21/components'/kind/(c['component_id']+'.tex')
        if hashlib.sha256(p.read_bytes()).hexdigest()!=c[field]:errors.append(f'{kind} identity: {c["component_id"]}')
assert len(checklist['components'])==215
assert len(checklist['c1_cells'])==45
assert len(checklist['c2_leaves'])==49

links=0
for p in HERE.glob('*.md'):
    for target in re.findall(r'\]\(([^\n)]+)\)', p.read_text(encoding='utf-8')):
        target=target.split('#',1)[0]
        if not target or re.match(r'^[a-z]+:',target):continue
        links+=1
        if not (p.parent/unquote(target)).exists():errors.append(f'broken link in {p.name}: {target}')

payload=json.loads((HERE/'SOURCE_MAP.json').read_text(encoding='utf-8'))
copied=0
for row in payload['records']:
    p=REPO/row['destination']
    if not p.exists() or hashlib.sha256(p.read_bytes()).hexdigest()!=row['sha256']:
        errors.append('source map mismatch: '+row['destination'])
    copied+=1

archive_results=[]
newdir=HERE/'materials/historical/analysis/exact-lift-2026-10-23/seven-extra-high-2026-10-26/fetched-for-publication'
for p in sorted(newdir.glob('*.zip')):
    with zipfile.ZipFile(p) as z:
        bad=z.testzip()
        if bad:errors.append(f'{p.name}: bad CRC for {bad}')
        archive_results.append({'file':p.name,'members':len(z.namelist()),'crc_ok':bad is None})

archive=HERE/'materials/historical/recovery/full_r21/archives/Q222_R21_Rederived_Proof_Package_2026-10-22.zip'
bindings=json.loads((HERE/'materials/historical/recovery/full_r21/records/PROOF_DEPENDENCIES.json').read_text(encoding='utf-8'))['bindings']
source_index=json.loads((HERE/'materials/historical/recovery/full_r21/records/COMPONENT_INDEX.json').read_text(encoding='utf-8'))
with zipfile.ZipFile(archive) as z:
    bad=z.testzip()
    if bad:errors.append('Original R21 archive CRC: '+bad)
    for row in bindings:
        path=source_index['inner_root']+'/'+row['path']
        try:data=z.read(path)
        except KeyError:
            errors.append('missing original archive binding: '+row['path']);continue
        if hashlib.sha256(data).hexdigest()!=row['sha256']:errors.append('archive binding hash: '+row['path'])

parsed=0
for p in newdir.glob('*.py'):
    try:ast.parse(p.read_text(encoding='utf-8-sig'),filename=p.name);parsed+=1
    except Exception as exc:errors.append(f'New checker syntax {p.name}: {exc}')

result={'scope':'delivery integrity only; scripts parsed, not executed; no new theorem certification',
        'verified_manifest_files':verified,'source_map_records_checked':copied,
        'component_source_and_presentation_hashes':430,'original_archive_dependency_bindings':len(bindings),
        'c1_cells':45,'c2_leaves':49,'local_links_checked':links,
        'new_original_archives':archive_results,'new_standalone_python_files_parsed':parsed,
        'errors':errors,'ok':not errors}
if '--write-receipt' in sys.argv:
    (HERE/'VERIFICATION.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print(json.dumps(result,indent=2))
raise SystemExit(bool(errors))

from pathlib import Path
import json,hashlib,shutil,subprocess,sys,os,time,ast
ROOT=Path(r'C:/Users/anubh/OneDrive/Documents/ChatGPT/ManuscriptWriting/q222/analytical-paper-2026-10-22')
SRC=ROOT/'github-q222-finite-sdp/independent-review/checks/gram-lift/ends'
OUT=Path(__file__).resolve().parent
import sympy
files=['check_c2_r0_face_k2.py','check_c2_n0_wall_k2.py','check_c1_zero_walls_k2.py','n0-wall/check_c2_n0_wall_family.py','n0-wall/data18.json']
r={'upstream_label':'3b1450d','python':sys.version,'sympy':sympy.__version__,'input_files':[],'runs':[],'environment':{}}
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','PYTHONDONTWRITEBYTECODE','PYTHONUTF8'):
    os.environ[key]='1';r['environment'][key]='1'
private=OUT/'private-copy'
private.mkdir(exist_ok=True)
for filename in files:
    p=SRC/filename;raw=p.read_bytes();dest=private/filename;dest.parent.mkdir(exist_ok=True,parents=True);shutil.copy2(p,dest)
    entry={'path':filename,'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()}
    if p.suffix=='.py':
        tree=ast.parse(raw)
        entry['assert_statements']=sum(isinstance(node,ast.Assert) for node in ast.walk(tree))
        entry['imports']=[ast.unparse(node) for node in tree.body if isinstance(node,(ast.Import,ast.ImportFrom))]
        for node in tree.body:
            if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='DATA' for t in node.targets):
                data=json.loads(ast.literal_eval(node.value.args[0])) if isinstance(node.value,ast.Call) else ast.literal_eval(node.value)
                entry['materialized_points']={k:len(v) for k,v in data.items()} if isinstance(data,dict) else len(data)
    else:
        data=json.loads(raw);entry['materialized_points']=len(data['pts']);entry['free_coordinates']=len(data['free'])
    r['input_files'].append(entry)
    print('INPUT',entry,flush=True)
def save(): (OUT/'results.json').write_text(json.dumps(r,indent=2),encoding='utf-8')
save()
for filename in files[:-1]:
    p=private/filename;log=OUT/(p.stem+'.log');start=time.perf_counter();print('START',filename,flush=True)
    with log.open('w',encoding='utf-8') as f:
        proc=subprocess.run([sys.executable,'-B','-u',p.name],cwd=p.parent,env=os.environ.copy(),stdout=f,stderr=subprocess.STDOUT,timeout=600)
    text=log.read_text();entry={'file':filename,'exit_code':proc.returncode,'elapsed_seconds':round(time.perf_counter()-start,6),'log':log.name,'pass_lines':sum('PASS (' in line for line in text.splitlines()),'final_line':text.splitlines()[-1] if text else ''}
    r['runs'].append(entry);save();print('END',entry,flush=True)
    if proc.returncode:print(text[-6000:],flush=True)
# Explicitly check no original was altered, including materialized data.
for entry in r['input_files']:
    entry['original_preserved']=hashlib.sha256((SRC/entry['path']).read_bytes()).hexdigest()==entry['sha256']
assert all(x['original_preserved'] for x in r['input_files'])
save()
print('BOUNDED ORIGINAL MATERIALIZED-WITNESS REPLAY COMPLETE',flush=True)


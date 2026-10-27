from pathlib import Path
import json,hashlib,shutil,subprocess,sys,os,time,ast
ROOT=Path(r'C:/Users/anubh/OneDrive/Documents/ChatGPT/ManuscriptWriting/q222/analytical-paper-2026-10-22')
REPO=ROOT/'github-q222-finite-sdp'
SRC=REPO/'independent-review/checks/gram-lift/ends/n0-opening'
OUT=Path(__file__).resolve().parent
import sympy
r={'upstream_label':'main2dc657f supplied by coordinator','python':sys.version,'sympy':sympy.__version__,'environment':{},'input_files':[]}
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','PYTHONDONTWRITEBYTECODE','PYTHONUTF8'):
    os.environ[key]='1';r['environment'][key]='1'
private=OUT/'private-copy'
if private.exists():raise RuntimeError('Private copy exists')
private.mkdir()
for name in ('verify_n0_opening.py','fixed_map.json','seed_wall_proposal.json','exact_results.json'):
    raw=(SRC/name).read_bytes();shutil.copy2(SRC/name,private/name)
    r['input_files'].append({'file':name,'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()})
proof=REPO/'independent-review/C2_N0_OPENING.md';raw=proof.read_bytes()
r['proof_identity']={'file':str(proof.relative_to(REPO)),'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()}
# Fix only the already delivered finite source; no dyadic/sign search.
orig=(private/'verify_n0_opening.py').read_text()
assert orig.count('for sign in (-1,1):')==1 and orig.count('for power in range(1,90):')==1
archived=json.loads((SRC/'exact_results.json').read_text())
assert archived['strict_step_t']=='1/16777216'
modified=orig.replace('for sign in (-1,1):','for sign in (1,):').replace('for power in range(1,90):','for power in (24,):')
(private/'verify_n0_opening.py').write_text(modified,encoding='utf-8')
r['private_execution_change']={'reason':'Remove adaptive source-step search; verify only the already materialized source +2^-24.','exact_replacements':[['for sign in (-1,1):','for sign in (1,):'],['for power in range(1,90):','for power in (24,):']],'executed_sha256':hashlib.sha256((private/'verify_n0_opening.py').read_bytes()).hexdigest(),'assert_statements':sum(isinstance(n,ast.Assert) for n in ast.walk(ast.parse(modified)))}
start=time.perf_counter()
with (OUT/'verify_n0_opening.log').open('w',encoding='utf-8') as f:
    p=subprocess.run([sys.executable,'-B','-u','verify_n0_opening.py'],cwd=private,env=os.environ.copy(),stdout=f,stderr=subprocess.STDOUT,timeout=300)
r['replay']={'exit_code':p.returncode,'elapsed_seconds':round(time.perf_counter()-start,6)}
print('REPLAY',r['replay'],flush=True)
if p.returncode:print((OUT/'verify_n0_opening.log').read_text(),flush=True)
assert p.returncode==0
r['output_comparisons']=[]
for name in ('exact_results.json','fixed_map.json'):
    before=(SRC/name).read_bytes();after=(private/name).read_bytes()
    compare={'file':name,'byte_identical':before==after,'newline_normalized_identical':before.replace(b'\r\n',b'\n')==after.replace(b'\r\n',b'\n'),'json_semantically_identical':json.loads(before)==json.loads(after)}
    r['output_comparisons'].append(compare);assert compare['json_semantically_identical'];print(compare,flush=True)
for entry in r['input_files']:
    entry['original_preserved']=hashlib.sha256((SRC/entry['file']).read_bytes()).hexdigest()==entry['sha256']
r['proof_identity']['original_preserved']=hashlib.sha256(proof.read_bytes()).hexdigest()==r['proof_identity']['sha256']
assert all(x['original_preserved'] for x in r['input_files']) and r['proof_identity']['original_preserved']
r['reported_pass_lines']=sum(line.startswith('PASS') for line in (OUT/'verify_n0_opening.log').read_text().splitlines())
r['derived_summary']={k:archived[k] for k in ('map_shape','rank','pivot_minor','jet_minor','unweighted_feature_count','weighted_feature_union_count','wall_theta_a_b_c_t_r','strict_theta_a_b_c_t_r','strict_weights_n_m_C_A_D_alpha','right_inverse_block_row_bounds','simple_coefficient_radius','strict_step_t')}
(OUT/'results.json').write_text(json.dumps(r,indent=2),encoding='utf-8')
print('PASS lines',r['reported_pass_lines'],'originals preserved',flush=True)
print(json.dumps(r['derived_summary'],indent=2),flush=True)

import ast,hashlib,json,os,pathlib,re,shutil,subprocess,sys,time,zipfile
root=pathlib.Path('C:/Users/anubh/OneDrive/Documents/ChatGPT/ManuscriptWriting/q222/analytical-paper-2026-10-22')
source=root/'analysis/exact-lift-2026-10-23/seven-extra-high-2026-10-26/fetched-1719'
package=source/'package'
work=pathlib.Path(__file__).parent
out=root/'analysis/c1_reflection_transport_exact_verification_2026_10_26_artifacts'
out.mkdir(exist_ok=True)
sha=lambda data:hashlib.sha256(data).hexdigest()
result={'archive':str(source/'original-reflection-package.zip'),'archive_sha256':sha((source/'original-reflection-package.zip').read_bytes()),'private_replay_dir':str(work),'python':sys.version,'integrity':{}}
before={p.name:sha(p.read_bytes()) for p in package.iterdir() if p.is_file()}
with zipfile.ZipFile(source/'original-reflection-package.zip') as archive:
    result['integrity']['zip_testzip_bad_member']=archive.testzip()
    result['integrity']['archive_members']=[{'name':z.filename,'size':z.file_size,'crc32':f'{z.CRC:08x}'} for z in archive.infolist()]
    rows=[]
    for z in archive.infolist():
        if z.is_dir():continue
        data=archive.read(z)
        path=package/z.filename
        rows.append({'name':z.filename,'archive_sha256':sha(data),'extracted_sha256':sha(path.read_bytes()) if path.is_file() else None,'exact_archive_extracted_match':path.is_file() and path.read_bytes()==data})
    result['integrity']['archive_vs_extracted']=rows
    sums_data=archive.read('SHA256SUMS.txt')
    rows=[]
    for line in sums_data.decode('utf8').splitlines():
        match=re.fullmatch(r'([0-9a-fA-F]{64})\s+\*?(.+)',line)
        if not match:raise ValueError(f'Invalid checksum line: {line!r}')
        expected,name=match.groups()
        archive_actual=sha(archive.read(name)); extracted_actual=sha((package/name).read_bytes())
        rows.append({'name':name,'expected':expected,'archive_actual':archive_actual,'extracted_actual':extracted_actual,'pass':expected==archive_actual==extracted_actual})
    result['integrity']['sha256sums_rows']=rows
    assert result['integrity']['zip_testzip_bad_member'] is None
    assert all(x['exact_archive_extracted_match'] for x in result['integrity']['archive_vs_extracted'])
    assert all(x['pass'] for x in rows)
    for name in ['verify_reflection_transport.py','transport_sparse_maps.json','EXACT_CHECKS.log']:
        (work/name).write_bytes(archive.read(name))
original_json=(work/'transport_sparse_maps.json').read_bytes()
script=(work/'verify_reflection_transport.py').read_text('utf8')
parsed=ast.parse(script)
result['checker_inspection']={'imports':[ast.unparse(n) for n in ast.walk(parsed) if isinstance(n,(ast.Import,ast.ImportFrom))],'source_sha256':sha(script.encode('utf8')),'assert_statement_count':sum(isinstance(n,ast.Assert) for n in ast.walk(parsed)),'literal_pass_print_count':sum(isinstance(n,ast.Constant) and isinstance(n.value,str) and n.value.startswith('PASS ') for n in ast.walk(parsed)),'file_output':'transport_sparse_maps.json written beside __file__','solver_calls_observed':False,'numerical_sweep_observed':False}
env=os.environ.copy()
for name in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS']:env[name]='1'
env['PYTHONHASHSEED']='0'
start=time.monotonic()
run=subprocess.run([sys.executable,str(work/'verify_reflection_transport.py')],cwd=work,env=env,capture_output=True,timeout=240)
result['replay']={'returncode':run.returncode,'elapsed_seconds':time.monotonic()-start,'stdout_sha256':sha(run.stdout),'stderr_sha256':sha(run.stderr),'stdout_exact_original_log_match':run.stdout==(work/'EXACT_CHECKS.log').read_bytes(),'stdout_normalized_original_log_match':run.stdout.replace(b'\r\n',b'\n')==(work/'EXACT_CHECKS.log').read_bytes().replace(b'\r\n',b'\n'),'generated_json_sha256':sha((work/'transport_sparse_maps.json').read_bytes()),'original_json_sha256':sha(original_json),'generated_json_byte_identical':(work/'transport_sparse_maps.json').read_bytes()==original_json,'generated_json_native_data_identical':json.loads((work/'transport_sparse_maps.json').read_text('utf8'))==json.loads(original_json)}
(work/'replay.stdout.log').write_bytes(run.stdout)
(work/'replay.stderr.log').write_bytes(run.stderr)
data=json.loads((work/'transport_sparse_maps.json').read_text('utf8'))
grid=data['extended_features']; base=data['oriented_Z16_features']; blocks=data['reflected_16_sparse_columns'];j0=data['Z6_to_ext_unmultiplied_indices'];j1=data['Z6_to_ext_times_v_indices']
assert len(grid)==364 and len(set(map(tuple,grid)))==364
assert grid==[[i,j,l] for i in range(7) for j in range(13) for l in range(4)]
assert len(base)==16 and len(blocks)==16 and sum(map(len,blocks))==32
assert all(isinstance(idx,int) and 0<=idx<364 and sign in (-1,1) for block in blocks for idx,sign in block)
assert len(j0)==len(j1)==196 and set(j0)|set(j1)==set(range(364))
assert all(grid[idx]==[i,j,l] for idx,(i,j,l) in zip(j0,[[i,j,l] for i in range(7) for j in range(7) for l in range(4)]))
assert all(grid[idx]==[i,j+6,l] for idx,(i,j,l) in zip(j1,[[i,j,l] for i in range(7) for j in range(7) for l in range(4)]))
result['native_json_checks']={'grid_rows':len(grid),'grid_unique_rows':len(set(map(tuple,grid))),'feature_columns':len(blocks),'signed_nonzero_entries':sum(map(len,blocks)),'j0_rows':len(j0),'j1_rows':len(j1),'embedding_union_coverage':len(set(j0)|set(j1)),'embedding_intersection':len(set(j0)&set(j1)),'coefficient_key_shape':data['coefficient_key_shape'],'checked_assertion_count_groups_metadata':data['checked_assertion_count_groups'],'pass_lines_observed':sum(line.startswith(b'PASS ') for line in run.stdout.splitlines()),'coefficient_linear_map_is_string_only':isinstance(data['coefficient_linear_map'],str),'all_checks_pass':True}
after={p.name:sha(p.read_bytes()) for p in package.iterdir() if p.is_file()}
result['integrity']['all_original_package_files_unchanged_after_replay']=before==after
assert before==after
for name in ['replay.stdout.log','replay.stderr.log','transport_sparse_maps.json']:
    shutil.copyfile(work/name,out/('replayed_transport_sparse_maps.json' if name=='transport_sparse_maps.json' else name))
(out/'integrity_and_replay_results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
print(json.dumps({'results_path':str(out/'integrity_and_replay_results.json'),'integrity_rows':len(result['integrity']['sha256sums_rows']),'replay':result['replay'],'native_json_checks':result['native_json_checks'],'original_files_unchanged':before==after},indent=2))
assert run.returncode==0
assert result['replay']['stdout_normalized_original_log_match']
assert result['replay']['generated_json_native_data_identical']

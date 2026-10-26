from pathlib import Path
import hashlib,json,zipfile
here=Path(__file__).resolve().parent
path=Path('C:/Users/anubh/Downloads/Q222_Global_Weighted_Source_and_Matrices.zip')
prefix='Q222_Global_Weighted_Certificate/'
with zipfile.ZipFile(path) as z:
    assert z.testzip() is None
    manifest=json.loads(z.read(prefix+'SHA256SUMS.json'))
    for name,sha in manifest.items():
        assert hashlib.sha256(z.read(prefix+name)).hexdigest()==sha,name
    reused=z.read(prefix+'reused/C1_UNIFORM_LOCAL_ORDER669.coo.json')
    local=here.parents[1]/'github-q222-finite-sdp/extracted-lifts/local-order669/root/C1_UNIFORM_LOCAL_ORDER669.coo.json'
    assert reused==local.read_bytes()
report={'archive':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
        'manifest_entries':len(manifest),'all_archive_hashes_match':True,
        'reused_order669_byte_identical_to_public':True,
        'scientific_status':'Fresh mathematical and matrix reviews pending; no full affine lift asserted.'}
(here/'pro133_integrity.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))

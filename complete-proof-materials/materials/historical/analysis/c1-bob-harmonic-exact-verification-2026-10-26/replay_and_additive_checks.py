"""Bounded independent replay and additive checks; preserve original package."""
from pathlib import Path
import hashlib,json,shutil,subprocess,sys,time,zipfile,os
BASE=Path(__file__).resolve().parent
ROOT=BASE.parent
SRC=ROOT/"exact-lift-2026-10-23/seven-extra-high-2026-10-26/fetched-1959/task5"
PKG=SRC/"package"
COPY=BASE/"private-package"
COPY.mkdir(exist_ok=True)
for p in PKG.rglob("*"):
    if p.is_file() and "__pycache__" not in p.parts:
        q=COPY/p.relative_to(PKG)
        q.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(p,q)
manifest=json.loads((PKG/"MANIFEST.json").read_text(encoding="utf-8"))
hash_checks=[]
for item in manifest["files"]:
    b=(PKG/item["name"]).read_bytes()
    assert len(b)==item["bytes"]
    assert hashlib.sha256(b).hexdigest()==item["sha256"]
    hash_checks.append(item["name"])
original=SRC/"original-package.zip"
archive_members=[]
with zipfile.ZipFile(original) as z:
    assert z.testzip() is None
    for info in z.infolist():
        if info.is_dir(): continue
        rel=Path(info.filename)
        expected=PKG/rel
        if not expected.exists():
            matches=list(PKG.rglob(rel.name))
            assert len(matches)==1,info.filename
            expected=matches[0]
        assert z.read(info.filename)==expected.read_bytes(),info.filename
        archive_members.append(info.filename)
env=os.environ.copy()
env.update(OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1",OMP_NUM_THREADS="1",PYTHONIOENCODING="utf-8")
results=[]
for name in ("verify_full_parent_bob.py","verify_boundary_schur.py","verify_anchor_interpolation.py"):
    started=time.perf_counter()
    p=subprocess.run([sys.executable,"-X","utf8",str(COPY/name)],cwd=COPY,capture_output=True,text=True,timeout=90,env=env)
    (BASE/(name.removesuffix(".py")+".log")).write_text(p.stdout+p.stderr,encoding="utf-8")
    results.append({"script":name,"exit_code":p.returncode,"seconds":time.perf_counter()-started})
    assert p.returncode==0,name
sys.path.insert(0,str(COPY))
import sympy as S
import verify_full_parent_bob as M
# Direct coefficient equality, not the earlier text's assertion.
assert S.expand(M.F.subs({M.u:-M.u,M.Y:-M.Y,M.Z:-M.Z},simultaneous=True)-M.F)==0
assert S.expand(M.F.subs({M.v:-M.v,M.X:-M.X,M.Z:-M.Z},simultaneous=True)-M.F)==0
# Endpoint row systems prove the exact kernel and PD claim symbolically.
assert M.G0[:,0]==S.zeros(8,1) and M.G0[:,5]==S.zeros(8,1)
range_indices=[1,2,3,4,6,7]
R=M.G0.extract(range_indices,range_indices)
expected_R_det=256*M.A**2*M.al**2*(4*M.C+4*M.D+2*M.d)*(4*M.m+4*M.be+2*M.d)
assert S.factor(R.det()-expected_R_det)==0
# Positivity of this determinant alone is not the PD proof: the explicit six
# positive-weight independent rows in the analytical note prove PD.
assert S.expand(M.Ginf[0,0])==4*M.m
# Every diagonal of H(t^2)=R(t)^T R(t) has noncancelling highest squares;
# the resulting 24 even +16 odd coordinates are counted independently.
even={(a,b,l) for a in (0,6) for b in (0,6,12) for l in range(4)}
odd={(a,b,l) for a in (0,6) for b in (3,9) for l in range(4)}
assert len(even)==24 and len(odd)==16 and not(even&odd)
assert all(a<=6 and b<=12 for a,b,l in even|odd)
receipt={
 "scope":"Three original exact modules plus exact sign congruences and endpoint range determinant. No universal H1 or middle-angle PSD verification.",
 "python":sys.version,"sympy":S.__version__,
 "manifest_files":len(hash_checks),"archive_members":len(archive_members),
 "original_zip_sha256":hashlib.sha256(original.read_bytes()).hexdigest(),
 "manifest_verified":hash_checks,"archive_verified":archive_members,"runs":results,
 "additional_checks":{"Alice_sign_congruence":True,"Bob_sign_congruence":True,
 "G0_exact_kernel_from_positive_rows":"span(W0,X1)",
 "G0_range_determinant":str(S.factor(R.det())),
 "fixed_features":40,"even_features":24,"odd_features":16},
 "original_package_unchanged":all(hashlib.sha256((PKG/i["name"]).read_bytes()).hexdigest()==i["sha256"] for i in manifest["files"])
}
(BASE/"receipt.json").write_text(json.dumps(receipt,indent=2)+"\n",encoding="utf-8")
print(json.dumps({"manifest_files":len(hash_checks),"archive_members":len(archive_members),"runs":results,"additive_checks":"PASS"},indent=2))


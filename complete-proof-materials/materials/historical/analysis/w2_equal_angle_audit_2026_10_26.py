"""Independent universal-word audit of the equal-angle ME C2 W2 obstruction."""
import json
from pathlib import Path
import subprocess
import sys
import sympy as s
from sympy.polys.matrices import DomainMatrix

r,z,c=s.symbols('r z c')
one=((),())
def red(a):
    st=[]
    for j in a:
        if st and st[-1]==j:st.pop()
        else:st.append(j)
    return tuple(st)
def word(a=(),b=()):return (tuple(a),tuple(b))
def mul(x,y):return (red(x[0]+y[0]),red(x[1]+y[1]))
def adj(x):return (tuple(reversed(x[0])),tuple(reversed(x[1])))
A0=word((0,));A1=word((1,));B0=word(b=(0,));B1=word(b=(1,))
ww=[one,A0,A1,B0,B1,mul(A0,A1),mul(A1,A0),mul(B0,B1),mul(B1,B0),
    mul(A0,B0),mul(A0,B1),mul(A1,B0),mul(A1,B1)]
h=[{A1:1,one:-1},{B1:1,one:-1},{mul(A0,B0):1,one:-1},
   {mul(A1,B1):1,one:-1},{B0:1,A0:-1},{mul(A0,B1):1,A0:-1},
   {mul(A1,B0):1,A0:-1},{mul(A0,A1):1,A0:-1},{mul(A1,A0):1,A0:-1},
   {mul(B0,B1):1,B0:-1},{mul(B1,B0):1,B0:-1}]
R=r*(1+z);p=(1+r)*(1-r*z);q=1+r*r*z;t=(r-1)*(1+r*z)
Kn=s.Matrix([[0,p,q,q,p,p,q],[0,-q,t,t,-q,-q,t],
 [-R,-p,-t,-q,q,q,-q],[R,0,0,0,0,0,0],[0,0,-R,-R,0,R,0],
 [0,R,0,0,0,0,0],[0,0,R,0,0,0,0],[0,0,0,R,0,0,0],
 [0,0,0,0,R,0,0],[0,0,0,0,0,R,0],[0,0,0,0,0,0,R]])
N=Kn/R
x=r*z;den0=r+x+2*r*x
src=[s.S.One,(1-r-x-3*r*x)/(4*den0),(1-r)*(1+x)/(4*den0),
 r*r*(1+x)/((1-r)*den0),(1+r)*(1+x)/(4*den0),
 (1+r)*(1-x)/(4*den0),x*x*(1+r)/((1-x)*den0)]
d,n,m,C,A,D,al=src
events=[(n,0,0,1,-1),(m,1,1,1,-1),(C,1,0,-1,1),
        (A,0,0,-1,1),(D,1,1,-1,1),(al,0,1,-1,-1)]
target={one:d,mul(A0,B0):-d/2,mul(A0,B1):d/2,mul(A1,B0):-d/2,mul(A1,B1):-d/2}
for weight,ix,iy,ax,by in events:
    aw=[A0,A1][ix];bw=[B0,B1][iy]
    for key,val in [(one,weight),(aw,ax*weight),(bw,by*weight),(mul(aw,bw),ax*by*weight)]:
        target[key]=target.get(key,0)+val

def pprod(a,b):
    out={}
    for wa,va in a.items():
        for wb,vb in b.items():
            w=mul(adj(wa),wb);out[w]=out.get(w,0)+va*vb
    return {w:s.expand(v) for w,v in out.items() if v!=0}
ks=[]
for j in range(7):
    out={}
    for i in range(11):
        for w,v in h[i].items():out[w]=out.get(w,0)+Kn[i,j]*v
    ks.append({w:s.expand(v) for w,v in out.items() if v!=0})
upper=list((i,j) for i in range(7) for j in range(i,7))
cols=[]
for i,j in upper:
    out=pprod(ks[i],ks[j])
    if i!=j:
        for w,v in pprod(ks[j],ks[i]).items():out[w]=s.expand(out.get(w,0)+v)
    cols.append(out)
words=sorted(set().union(*(set(t) for t in cols),set(target)),key=lambda w:(len(w[0])+len(w[1]),w))
AA=s.Matrix([[cc.get(w,0) for cc in cols] for w in words])
SDen=4*r*(1-r)*(1-r*z)*(1+z+2*r*z)
rhs=s.Matrix([s.cancel(R*R*SDen*target.get(w,0)) for w in words])
assert all(s.denom(v)==1 for v in rhs)

if '--rref' in sys.argv or '--sparse' in sys.argv:
    trace=[]
    if '--sparse' in sys.argv:
        free_wanted={upper.index(ij) for ij in [(3,3),(4,5),(4,6),(5,5),(6,6)]}
        pending=[(i,{j:v for j,v in enumerate(list(AA.row(i))+[rhs[i]]) if v!=0}) for i in range(len(words))]
        solved=[]; pivots=[]
        allowed=[r,1-r,1-r*z,1+r,1+r*z,1+z,1+z+2*r*z]
        while True:
            choices=[(len(row),s.count_ops(row[j]),len(str(row[j])),ii,j)
                     for ii,(label,row) in enumerate(pending)
                     for j in row if j<28 and j not in free_wanted]
            if not choices:break
            _,_,_,ii,j=min(choices)
            label,row=pending.pop(ii); pivot=s.factor(row[j])
            trace.append({'word':str(words[label]),'entry':[upper[j][0]+1,upper[j][1]+1],'pivot':str(pivot)})
            prow={k:s.cancel(v/pivot) for k,v in row.items()}
            for ii,(lab,old) in enumerate(pending):
                if j not in old:continue
                cf=old[j];new={}
                for k in old.keys()|prow.keys():
                    if k==j:continue
                    vv=s.cancel(old.get(k,0)-cf*prow.get(k,0))
                    if vv!=0:new[k]=vv
                pending[ii]=(lab,new)
            solved.append(prow);pivots.append(j)
        assert all(not row for label,row in pending),str(pending)
        rr=s.zeros(len(pivots),29);common=1
        # Back substitution retains the chosen five free coordinates.
        for row in range(len(solved)-1,-1,-1):
            pr=solved[row]
            for later in range(row+1,len(solved)):
                jj=pivots[later]
                if jj in pr:
                    cf=pr.pop(jj)
                    for k,v in solved[later].items():
                        if k==jj:continue
                        pr[k]=s.cancel(pr.get(k,0)-cf*v)
            solved[row]={k:v for k,v in pr.items() if v!=0}
            for k,v in solved[row].items():rr[row,k]=v
    else:
        DM=DomainMatrix.from_Matrix(AA.row_join(rhs)).convert_to(s.ZZ.poly_ring(r,z))
        rr,common,pivots=DM.rref_den(method='FF')
        rr=rr.to_Matrix();common=DM.domain.to_sympy(common)
    free=[i for i in range(28) if i not in pivots]
    vv=s.symbols('free0:'+str(len(free)))
    sol={i:vv[j] for j,i in enumerate(free)}
    for row,pi in enumerate(pivots):
        assert pi<28
        sol[pi]=s.cancel(rr[row,28]/common-sum(rr[row,j]*sol[j]/common for j in free))
    S=s.zeros(7)
    for k,(i,j) in enumerate(upper):S[i,j]=S[j,i]=sol[k]/SDen
    req=[S[0,0]-S[4,5]+2*S[4,6]+S[6,6],S[4,4]+S[4,5],
         S[0,4]+S[0,6],S[0,3],S[3,4],S[0,5],S[5,6],
         S[3,6]+S[6,6],
         S[2,4]-r*(1+r)*(1+z)**2/(32*(1-r)*(1+z+2*r*z))]
    assert all(s.cancel(v)==0 for v in req)
    # These polynomials do not vanish on the stated real domain.
    allowed=[r,1-r,1-r*z,1+r,1+r*z,1+z,1+z+2*r*z]
    factors=set();pivot_factors=set()
    for item in trace:
        pivot=s.sympify(item['pivot'])
        for part in s.fraction(pivot):
            for factor,exp in s.factor_list(part)[1]:
                pp=s.Poly(factor,r,z)
                domain_nonzero=(any(pp.monic()==s.Poly(a,r,z).monic() for a in allowed)
                    or all(t>0 for t in pp.coeffs()) or all(t<0 for t in pp.coeffs()))
                assert domain_nonzero, 'Uncertified pivot locus: '+str(factor)
                pivot_factors.add(str(factor))
    for v in sol.values():
        denominator=s.denom(s.cancel(v/SDen))
        for factor,exp in s.factor_list(denominator)[1]:
            assert any(s.Poly(factor,r,z).monic()==s.Poly(a,r,z).monic() for a in allowed),str(factor)
            factors.add(str(factor))
    print(json.dumps({'word_rows':len(words),'coefficient_rank':len(pivots),
        'free_entries_one_based':[[upper[j][0]+1,upper[j][1]+1] for j in free],
        'all_nine_essential_relations':True,'denominator_factors':sorted(factors),
        'all_pivots_nonzero_throughout_domain':True,'pivot_factors':sorted(pivot_factors),
        'triangular_elimination':trace,
        'S35':str(r*(1+r)*(1+z)**2/(32*(1-r)*(1+z+2*r*z)))}))
    raise SystemExit()

checks={}
# Thirteen independent group words, and the two-evaluation kernel is exactly h.
ev=lambda w,a:s.prod(a[i] for i in w[0])*s.prod(a[2+i] for i in w[1])
evals=s.Matrix([[ev(w,a) for w in ww] for a in [(1,1,1,1),(-1,1,-1,1)]])
rootcols=s.Matrix([[hi.get(w,0) for hi in h] for w in ww])
checks['full_W2_two_classical_zero_compression']=len(set(ww))==13 and evals.rank()==2 and rootcols.rank()==11 and evals*rootcols==s.zeros(2,11)
Zp=s.diag(1,-1);Xp=s.Matrix([[0,1],[1,0]]);Id=s.eye(2)
Mc=((1-c*c)*Zp+2*c*Xp)/(1+c*c)
opsA=[Zp,Mc];opsB=[Zp,Mc]
def mat(w):
    a=Id;b=Id
    for k in w[0]:a=a*opsA[k]
    for k in w[1]:b=b*opsB[k]
    return s.kronecker_product(a,b)
psi=s.Matrix([1,-r*c,r*c,1])
T=s.Matrix.hstack(*[sum((v*mat(w) for w,v in hi.items()),s.zeros(4))*psi for hi in h])
checks['printed_kernel_all_28_entries']=all(s.cancel(v)==0 for v in T*N.subs(z,c*c))
checks['rank_four_contact_pivot']=s.factor(T[:,[0,1,2,4]].det(method='domain-ge'))==-64*r**3*c**4/(1+c*c)
checks['seven_identity_free_rows']=N[[3,5,6,7,8,9,10],:]==s.eye(7)
# Recover the equal-angle coefficients directly from the original native formulas.
Ac=1+z*r;tt=1-r;hh=r;gg=jj=R
Delta=(2*z+4*z*r+2*(1+z)**2) # replace below with the unsimplified original formula
Delta=(1+z)*z+2*z*(1+z)*r+(1+z)*(1+z)+2*z*(1+z)*r+1+z
NN=(1+z)*(1+r)*(2-tt)-(1+z)*(1-r)*tt
native=[1,((1+r)*Ac*(Ac+z*tt)-r*Delta)/(2*r*Delta),
 (1+z)*(1-r)*Ac*(Ac+z*tt)/(2*gg*Delta),Ac*NN/(2*tt*Delta),
 Ac*(1+r)*(Ac+z*tt)/(2*hh*Delta),
 (1+z)*(1+r)*(1-z*r)*(Ac+z*tt)/(2*jj*Delta),
 z*z*(1+r)*NN/(2*(1-z*hh)*Delta)]
checks['surface_source_equals_original_native_formulas']=all(s.cancel(v-w)==0 for v,w in zip(src,native))
checks['strict_seed_source']= [s.factor(v.subs({r:s.Rational(1,4),z:1})) for v in src]==[1,s.Rational(1,8),s.Rational(3,8),s.Rational(1,6),s.Rational(5,8),s.Rational(3,8),s.Rational(1,6)]
checks['n_zero_seed_source']=[s.factor(v.subs({r:s.Rational(1,3),z:1})) for v in src]==[1,0,s.Rational(1,4),s.Rational(1,4),s.Rational(1,2),s.Rational(1,4),s.Rational(1,4)]
S35=r*(1+r)*(1+z)**2/(32*(1-r)*(1+z+2*r*z))
checks['both_S35_seed_values']=[S35.subs({r:s.Rational(1,4),z:1}),S35.subs({r:s.Rational(1,3),z:1})]==[s.Rational(1,48),s.Rational(1,32)]
checks['n_boundary_derivative']=s.diff(src[1],r).subs({r:s.Rational(1,3),z:1})==-s.Rational(9,8)
checks['unitary_distinct_words_trace_bound']=all((mul(adj(wi),wj)==one)==(i==j) for i,wi in enumerate(ww) for j,wj in enumerate(ww))
worker=subprocess.run([sys.executable,__file__,'--sparse'],capture_output=True,text=True,timeout=30)
if worker.returncode:
    raise RuntimeError(worker.stderr+worker.stdout)
rr=json.loads(worker.stdout)
checks['symbolic_universal_elimination']=rr['all_nine_essential_relations'] and rr['word_rows']==41 and rr['coefficient_rank']==23
result={'checks':checks,'all_pass':all(checks.values()),'universal_elimination':rr,
 'scope':'Full 13-word W2 exclusion on r,x>0,r+x+3rx<=1; x=r*z. Unequal-angle positive Gram is not audited here.',
 'worker_timeout_seconds':30,'no_solver_or_parameter_sweep':True}
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print(json.dumps(result,indent=2))
if not result['all_pass']:raise SystemExit(1)

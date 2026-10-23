import sympy as S,pickle,json
from pathlib import Path
from fractions import Fraction as F
O=Path(__file__).resolve().parents[2]/'data'
Q=pickle.load((O/'lII_generic.pkl').open('rb'))['Q']
u,v,z,k,H,B,r,cp,g,q,ss,ls,als,al,nu,mu,aa=S.symbols('u v z k H B r cp g q ss ls als al nu mu aa')
hh=S.symbols('hh',positive=True)
tt,zz,kk,db,ds,dc,dg,dal,dA,dm,dn,dq=S.symbols('tt zz kk db ds dc dg dal dA dm dn dq')
ff=S.Matrix(S.symbols('Nn Xx Zz Hh Rr Tr'));Nn,Xx,Zz,Hh,Rr,Tr=ff
con=Tr**2-u*u*Rr**2-(u*u-z)/ss*Hh*Tr+u*(u*u-z)/ss*Rr*(u*Xx-Zz)
Qc=S.hessian(con,list(ff))/2
c0=1/(2*(1+hh*hh))
sub={H:hh*hh,B:hh*hh*(1+db),ss:(1+ds)/(hh*hh),cp:hh*hh*(c0+dc),g:hh*hh*(c0+dg),q:tt+dq,ls:hh*hh/4,als:hh**4/8,al:hh**4*(c0+dal),nu:c0-r/2+dn,mu:((1-r)*c0+dm)/hh**4,aa:(c0+dA)/(hh*hh),z:hh*hh*zz,k:hh*kk}
Q=Q+4/H*Qc
T=S.diag(1/hh,1,1/(hh*hh),1/hh,hh,hh)
T[1,2]=-u/(hh*hh);T[1,4]=-1/(hh*(1+hh*hh))
Q=T.T*Q.subs(sub,simultaneous=True)*T
zero={x:0 for x in [zz,kk,db,ds,dc,dg,dal,dA,dm,dn,dq]}
base=Q.subs(zero).subs({u:0,v:0})
variables=[u,v,hh,r,tt,zz,kk,db,ds,dc,dg,dal,dA,dm,dn,dq]
eps=F(1,2**64)
bounds={u:F(1,2**12),v:F(1,2**25),r:F(1),tt:F(2),zz:eps,kk:F(1,2**32),db:2*eps,ds:4*eps,dc:16*eps,dg:64*eps,dal:16*eps,dA:16*eps,dm:16*eps,dn:16*eps,dq:64*eps}
hmin=F(1,64)
def rational_bound(expr):
    num,den=S.fraction(S.cancel(expr))
    den=S.factor(den)
    terms=S.factor_list(den)[1]
    const=F(S.factor_list(den)[0]);hp=0;denlower=const
    for p,e in terms:
        if p==hh:hp+=e
        elif p==hh**2+1:pass
        elif p in [db+1,ds+1]:denlower*=F(1,2)**e
        else:raise ValueError(('unexpected denom',den,p,e))
    P=S.Poly(num,*variables)
    total=F(0)
    for mon,co in P.terms():
        term=abs(F(co))
        for var,power in zip(variables,mon):
            if var==hh:
                power-=hp
                term*=hmin**power if power<0 else F(1)
            else:term*=bounds[var]**power
        total+=term
    return total/denlower,len(P.terms())
records=[];source_sum=F(0);testing_sum=F(0)
for i in range(6):
 for j in range(i,6):
    print('entry',i,j,flush=True)
    entry=S.cancel(Q[i,j])
    lead=S.cancel(entry.subs(zero))
    source=S.cancel(entry-lead)
    testing=S.cancel(lead-base[i,j])
    sb,ns=rational_bound(source);tb,nt=rational_bound(testing)
    mult=1 if i==j else 2
    source_sum+=mult*sb;testing_sum+=mult*tb
    records.append({'i':i,'j':j,'source_bound':str(sb),'testing_bound':str(tb),'source_terms':ns,'testing_terms':nt})
    Q[i,j]=Q[j,i]=entry
print('source sum',float(source_sum),'testing sum',float(testing_sum))
receipt={'source_sum':str(source_sum),'testing_sum':str(testing_sum),'sum_lt_one_quarter':source_sum+testing_sum<F(1,4),'records':records}
(O/'lII_uniform_bound.json').write_text(json.dumps(receipt,indent=2))
with (O/'lII_normalized_full.pkl').open('wb') as fp:pickle.dump({'Q':Q,'base':base,'variables':variables,'bounds':bounds,'source_substitution':sub},fp)

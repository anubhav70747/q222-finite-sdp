"""Independent rational TL3 reconstruction; no sampling/solver."""
import json
from pathlib import Path
import sympy as S

def main():
    a,B,h,t,r=S.symbols('a B h t r',nonzero=True)
    e=1-r;s=1/a;x=a*B;y=t*B;l=1-h;K=1+x;ca=1+x*r;Sc=1+B*l+x*r
    Qc=1-h*B;E=1-t*x*r;R=1+x*h;g=R-a*e;j=ca-l/t
    De=a+1+x+(1+x)*B+2*x*(1+x)*r+a**3*B**2*r**2+2*a*a*B*r+t*B*(1+a+x+a*a*B*r*(2-r))
    M=K*(1+a*r)*(h+t)-a*(1+y)*e*l
    N=(1+a*r)*ca*Sc-r*De
    # Original native C2 formulas, using squared a0 to avoid radicals.
    aa=t/a;z=x;bb=1/a;tt=l/t
    hh=bb-aa*tt;gg=(1+z)*bb-(1-r+aa*z*tt);jj=1+z*r-tt
    A0=1+z*r
    delta=(1+z)*z*bb**2+2*z*(1+z)*bb*r+(1+z)*(1+aa*z)*bb+z*z*(1-aa)*r*r+2*z*(1+aa*z)*r+1+aa*z
    NN=(1+z)*(bb+r)*(bb+aa*(1-tt))-aa*(1+aa*z)*(1-r)*tt
    src2=[((bb+r)*A0*(A0+aa*z*tt)-r*delta)/(2*r*delta),
          (1+aa*z)*(1-r)*A0*(A0+aa*z*tt)/(2*gg*delta),A0*NN/(2*aa*tt*delta),
          aa*A0*(bb+r)*(A0+aa*z*tt)/(2*hh*delta),
          (1+z)*(bb+r)*(bb-aa*z*r)*(A0+aa*z*tt)/(2*aa*jj*delta),z*z*(bb+r)*NN/(2*(1-z*hh)*delta)]
    short2=[N/(2*r*De),a*a*(1+y)*e*ca*Sc/(2*g*De),ca*M/(2*l*De),
            t*ca*(1+a*r)*Sc/(2*h*De),K*(1+a*r)*E*Sc/(2*t*j*De),B*B*(1+a*r)*M/(2*Qc*De)]
    zero=lambda q:S.cancel(q)==0
    checks={}
    for i,name in enumerate(['n','m','C','A','D','alpha']):checks['native_C2_'+name]=zero(src2[i]-short2[i])
    checks['native_denominator']=zero(a*delta-De)
    checks['native_numerator']=zero(a*a*NN-M)
    # Native C1 NC2/NC4, not a physicality assertion.
    eta=t;zeta=s;r1=h;theta=e;h1=r;a1=l/t;t1=a*e;Y=x*eta*zeta
    RR=1+x*r1;TT=1+x*h1;GG=h1+zeta;MM=RR-t1;JJ=TT-a1;EE=1-x*eta*h1;QQ=1-x*zeta*r1
    P4=TT*EE+x*(1-r1)*GG
    Pa=RR*h1*(1+Y)+zeta*MM*EE
    Pd=1+2*x*h1+r1*x*x*GG+x*zeta+eta*x*(K*zeta+x*h1*theta)
    src1=[K*x*(r1+eta)*h1*(1+Y)/(2*MM*zeta*Pd),K*(r1+eta)*EE/(2*(1-r1)*Pd),
          eta*P4/(2*r1*Pd),(K*(1+Y)*(r1+eta)-(eta*TT-1+r1)*Pd)/(2*JJ*eta*Pd),
          x*x*zeta*(r1+eta)*Pa/(2*QQ*Pd),(K*x*(r1+eta)*GG-theta*Pd)/(2*theta*Pd)]
    pi=K*B*(h+t)*(1+a*r)/(e*Pd);kap=(1+2*short2[0])/pi;Th=short2[0]/B
    Pk=ca*Sc+K*B*(h+t);PD=ca*t*j-K*(h+t);Pal=(r+s)*ca-e*Qc
    for name,expr in [('bridge2',e*Pd-K*B*(h+t)*(1+a*r)-N),('bridge3',P4-r*Pd-N),('bridge4',e*Sc*E-r*B*M-N),
                      ('short_Pd',Pd-ca*(1+y)-B*(1+a*r)*(K*E+x*t*j)),
                      ('short_Palpha',ca*Pa-E*Pal-r*B*M),('beta_slack',src1[5]+(1-pi)/2),
                      ('kappa',kap-1-2*Th*Pk/(K*(h+t)))]:checks[name]=zero(expr)
    expected=[0,-Th*ca/l,-Th*t*ca*Sc/(K*(h+t)),Th*Sc*PD/(K*(h+t)*t*j),-short2[0]*x*Pal/(K*Qc),1-kap,-short2[0]]
    actual=[short2[1]-kap*src1[0],short2[2]-kap*src1[1],h*(short2[3]-kap*src1[2]),
            short2[4]-kap*src1[3],short2[5]-kap*src1[4],1-kap,-kap*src1[5]+(1-kap)/2]
    for name,q,w in zip(['m','C','Ah','D','alpha','d','paired_beta'],actual,expected):checks['TL3_'+name]=zero(q-w)
    # Literal R10.2, keeping b,w unrestricted until the exact comparison.
    ss,ww,bbig=S.symbols('ss ww bbig',nonzero=True)
    sub={a:1/ss,B:ss**2*bbig,r:1-ss*ww}
    xx=ss*bbig;rr=1-ss*ww
    Qstar=1+ss*xx*(1+t)+xx*(h+2*rr+t)+ss*xx**2*(h+t)+xx**2*(h*(rr+1)+t*(1+rr-rr**2))
    Kstar=(1+ss)*(h+t)*(1+xx)
    Dstar=1+ss*(1+xx)+ss**2*xx*(1+xx)+2*ss*xx*(1+xx)*rr+xx**2*rr**2+2*xx*rr+t*(ss**2*xx*(1+xx)+ss*xx+ss*xx**2*rr*(2-rr))
    checks['R10_exact_numerator']=zero(N.subs(sub,simultaneous=True)-ss*(ww*Qstar-bbig*Kstar))
    checks['R10_exact_denominator']=zero(ss*De.subs(sub,simultaneous=True)-Dstar)
    ee,dd=S.symbols('E delta',positive=True)
    CA=S.Matrix([[4*(ee*ee+4)/ee,4*(ee+2)/(ee*dd)],[4*(ee+2)/(ee*dd),4*(ee+4)/(ee*ee*dd*dd)]])
    SE=S.Matrix([[8*(ee*ee+2),8*(ee+1)],[8*(ee+1),8*(ee+2)/ee]])
    CT=S.Matrix([[1,-S.Rational(1,2),0],[-S.Rational(1,2),1,S.Rational(1,2)],[0,S.Rational(1,2),1]])
    checks['CA_residue_determinant']=zero(CA.det()-256/(ee**3*dd**2))
    checks['SE_residue_determinant']=zero(SE.det()-64*(1+4/ee))
    checks['circuit_leading_minors']=[CT[:i,:i].det() for i in range(1,4)]==[1,S.Rational(3,4),S.Rational(1,2)]
    assert all(checks.values())
    out={'checks':checks,'all_pass':True,'scope':'Exact source/TL3 identities and displayed strict residue matrices only; not a replay of the CA6/SE9 scalar coefficient lift.'}
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()

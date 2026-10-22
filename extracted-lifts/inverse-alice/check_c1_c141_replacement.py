"""Exact polynomial identities for the complete C14.1 three-page proof."""
import json
from pathlib import Path
import sympy as s

out={}
def zero(name,expr):
    value=s.factor(s.together(expr))
    assert value==0,(name,value)
    out[name]=True

p,w,W,X,Y,Z,q,b,f=s.symbols('p w W X Y Z q b f')
m,C,A,D,alpha,beta=s.symbols('m C A D alpha beta')
z=q*q;Q=1-b;r=b*z;H0=1+z*Q;B=H0/(b*Q)
C0=z/(2*(1-r));A0=H0/(2*b*z);al0=1/(2*Q*z)
Dp=1+p*p;Dw=1+w*w
E=q*W;F=Y/q;G=w*Z/q;R=X+w*W;S=p*X;J=F-b*G
K=R*R+Dw*S*S+z*(G-F)**2+p*p*Dw*z*F*F-2*p*E*(G+w*w*F)+2*q*w*S*(G-F)
eventM=p*W-p*w*X+Y-w*Z
eventD=w*W+X-p*w*Y-p*Z
raw=(4*m*eventM**2+4*C*Dw*(W-p*Y)**2+4*A*Dp*Dw*Y**2
     +(4*D+2)*eventD**2+(4*alpha+2)*Dp*(Y-w*Z)**2
     +(4*beta+2)*Dw*(p*X+Z)**2-2*Dp*Dw*Z**2)
events=(4*m*eventM**2+4*C*Dw*(W-p*Y)**2+4*A*Dp*Dw*Y**2
        +4*D*eventD**2+4*alpha*Dp*(Y-w*Z)**2+4*beta*Dw*(p*X+Z)**2)
zero('IA2_full_core',raw-events-2*K)
leadingLeft=(2*z*(G-F)**2+2*p*p*z*F*F-4*p*E*G
             +2/(1-r)*(E-p*z*F)**2+2*H0/b*Dp*F*F+2/Q*Dp*(F-G)**2)
leadingRight=2/(1-r)*(E-p*z*J-p*G)**2+2*B*Dp*J*J+2*H0*G*G
zero('IA3_exact_three_square_cancellation',leadingLeft-leadingRight)
zero('IA4_contraction',1-z*z/B-1/H0-z*Q*(1-r)/H0)

excess=(4*m*eventM**2+4*D*eventD**2+4*beta*Dw*(p*X+Z)**2
        +4*(C-f*C0)*Dw*(W-p*Y)**2+4*(A-f*A0)*Dp*Dw*Y**2
        +4*(alpha-f*al0)*Dp*(Y-w*Z)**2)
Pi=(excess+2*(1-f)*(R*R+S*S+z*(G-F)**2+p*p*z*F*F)
    +2*w*w*S*S+2*p*p*w*w*z*F*F+4*f*C0*w*w*(W-p*Y)**2
    +4*f*A0*Dp*w*w*Y**2)
Qlead=2*R*R+2*S*S+leadingRight
error=-4*(1-f)*p*E*G-4*p*w*w*E*F+4*q*w*S*(G-F)
zero('IA5_IA6_complete_parent',raw-f*Qlead-Pi-error)

P,t,Es,Js,Gs,Rg,Sg,a,c=s.symbols('P t Es Js Gs Rg Sg a c')
L=a*Js+c*Gs
Qnorm=2*Rg**2+2*Sg**2+2*P**2/(1-r)*(Es-t*L)**2+2*(1+P**2*t*t)*Js**2+2*Gs**2
Bnorm=2*Rg**2+2*Sg**2+2*(Es-t*L)**2+2*(1+t*t)*Js**2+2*Gs**2
zero('IA7_exact_clipping',Qnorm-Bnorm
     -2*(P**2/(1-r)-1)*(Es-t*L)**2-2*(P**2-1)*t*t*Js**2)
reserve=(s.Rational(3,2)*(Rg*Rg+Sg*Sg)+s.Rational(3,2)*(Es-4*t*L/3)**2
         +s.Rational(3,2)*Js**2+s.Rational(5,6)*Gs**2
         +s.Rational(2,3)*(1-t*t)*Gs**2+s.Rational(4,3)*t*t*Js**2
         +s.Rational(2,3)*t*t*((1-a*a-c*c)*(Js*Js+Gs*Gs)+(c*Js-a*Gs)**2))
zero('IA8_exact_reserve',Bnorm-(Rg*Rg+Sg*Sg+Es*Es+Js*Js+Gs*Gs)/2-reserve)

invB,invH=s.symbols('invB invH')
En=P*Es;Fn=invB*Js+b*invH*Gs;Gn=invH*Gs;pn=P*t
errnorm=-4*(1-f)*pn*En*Gn-4*pn*w*w*En*Fn+4*q*w*Sg*(Gn-Fn)
roots=s.Matrix([Rg,Sg,Es,Js,Gs]);Emat=s.zeros(5)
Emat[2,4]=Emat[4,2]=-2*(1-f)*P*P*t*invH-2*P*P*b*t*w*w*invH
Emat[2,3]=Emat[3,2]=-2*P*P*t*w*w*invB
Emat[1,3]=Emat[3,1]=-2*q*w*invB
Emat[1,4]=Emat[4,1]=2*q*Q*w*invH
zero('IA9_all_normal_entries',errnorm-(roots.T*Emat*roots)[0])

path=Path(__file__).resolve().with_name('c1_c141_replacement_checks.json')
path.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n',encoding='utf-8')
print(json.dumps(out,indent=2))

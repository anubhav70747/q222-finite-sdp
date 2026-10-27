#!/usr/bin/env python3
"""Exact source-only C1 contact/derivative skew recipe.
Input: exact symmetric native P0,P1,P2, source contact uc,vc,psi, and
contact scalar tau chosen from the strict source-only three-minor interval.
Output: quadratic polynomial cross-skew T(v), native contact Schur margin.
All inversions are pure source matrices with proved positive determinants.
No SDP, no testing-dependent coefficient, no global PSD claim.
"""
import sympy as s

def skew4(T):
    S=s.zeros(4)
    S[:2,2:]=T
    S[2:,:2]=-T.T
    return S

def skew8(T):
    S=skew4(T); J=s.zeros(8)
    J[:4,4:]=S; J[4:,:4]=-S
    return J

def H_matrix(v,P0,P1,P2,Tpoly):
    Dv=1+v**2; S=skew4(Tpoly)
    return (Dv*P0).row_join(Dv*P1/2+S).col_join((Dv*P1/2-S).row_join(Dv*P2))

def source_contact_jet(v,P0,P1,P2,uc,vc,psi,tau,verify=True):
    """Symbolic recipe, with supplied exact/algebraic native source data only.
    Requires the actual C1 contact six-jet and a tau making both contact
    blocks strictly PD on the 3-dimensional contact complement.
    """
    w=s.Matrix(psi[:2]);y=s.Matrix(psi[2:]);sw=w.dot(w);sy=y.dot(y)
    Dstar=1+vc**2
    A0=(Dstar*P0[:2,:2]).subs(v,vc)
    C2=(Dstar*P2[2:,2:]).subs(v,vc)
    L=(Dstar*P1[:2,2:]/2).subs(v,vc)
    p= -A0*w/uc-L*y
    q= -L.T*w-uc*C2*y
    compat=s.cancel(w.dot(p)-y.dot(q))
    if compat!=0: raise ValueError('native contact parity compatibility fails: '+str(compat))
    wp=s.Matrix([-w[1],w[0]]);yp=s.Matrix([-y[1],y[0]])
    Tstar=p*y.T/sy+w*q.T/sw-w.dot(p)*w*y.T/(sw*sy)+tau*wp*yp.T
    Tstar=Tstar.applyfunc(s.cancel)
    Hfixed=H_matrix(v,P0,P1,P2,Tstar)
    M=Hfixed.subs(v,vc)
    xA=s.Matrix([*w,0,0,0,0,*(uc*y)])
    xB=s.Matrix([0,0,*y,*(uc*w),0,0]);xp=xA+xB
    if verify:
        assert (M*xA).applyfunc(s.cancel)==s.zeros(8,1)
        assert (M*xB).applyfunc(s.cancel)==s.zeros(8,1)
    PiN=xA*xA.T/xA.dot(xA)+xB*xB.T/xB.dot(xB)
    Msharp=(M+PiN).inv()-PiN
    Msharp=Msharp.applyfunc(s.cancel)
    if verify:
        assert (M*Msharp+s.zeros(8)- (s.eye(8)-PiN)).applyfunc(s.cancel)==s.zeros(8)
    H1base=Hfixed.diff(v).subs(v,vc)
    Tunit=w*y.T/(sw*sy)
    coeff=s.cancel(-(xA.T*H1base*xA)[0]/(2*uc))
    Tfirst=coeff*Tunit
    H1=H1base+skew8(Tfirst)
    if verify:
        assert s.cancel((xA.T*H1*xA)[0])==0
        assert s.cancel((xB.T*H1*xB)[0])==0
    g1=skew8(wp*y.T)*xp
    g2=skew8(w*yp.T)*xp
    R=s.Matrix([[s.cancel((x.T*Msharp*y)[0]) for y in [g1,g2]] for x in [g1,g2]])
    b=H1*xp
    h=s.Matrix([s.cancel((x.T*Msharp*b)[0]) for x in [g1,g2]])
    optimal=(-R.inv()*h).applyfunc(s.cancel)
    Tfirst=Tfirst+optimal[0]*wp*y.T+optimal[1]*w*yp.T
    Tfirst=Tfirst.applyfunc(s.cancel)
    H1=H1base+skew8(Tfirst)
    bopt=H1*xp
    H2base=Hfixed.diff(v,2).subs(v,vc)/2
    H2base=H2base.applyfunc(s.cancel)
    mixed=H2base-H1*Msharp*H1
    kappa=s.cancel((xp.T*mixed*xp)[0])
    kA=s.cancel((xA.T*mixed*xA)[0])
    shift=s.cancel((kappa/2-kA)/(2*uc))
    Tsecond=(shift*Tunit).applyfunc(s.cancel)
    H2=H2base+skew8(Tsecond)
    if verify:
        assert s.cancel((xA.T*(H2-H1*Msharp*H1)*xA)[0]-kappa/2)==0
        assert s.cancel((xB.T*(H2-H1*Msharp*H1)*xB)[0]-kappa/2)==0
        assert s.cancel((bopt.T*Msharp*g1)[0])==0
        assert s.cancel((bopt.T*Msharp*g2)[0])==0
    Tpol=(Tstar+(v-vc)*Tfirst+(v-vc)**2*Tsecond).applyfunc(s.cancel)
    Hpol=H_matrix(v,P0,P1,P2,Tpol)
    return dict(Tpol=Tpol,Hpol=Hpol,Hcontact=M,Msharp=Msharp,kerA=xA,kerB=xB,
                contact_skew=Tstar,first_skew=Tfirst,second_skew=Tsecond,
                physical_schur=kappa,gram_2x2=R)

if __name__=='__main__':
    print('source_contact_jet is a symbolic function. Supply native source-only inputs and a strict algebraic contact-gauge tau.')
    print('The universal algebraic identities are tested in verify_contact_geometry.py and verify_two_schur.py.')

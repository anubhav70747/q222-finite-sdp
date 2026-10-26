"""Exact rational native C2 witness and universal-word obstruction.

Self-contained; requires SymPy. No SDP, general determinant expansion,
parameter scan, browser result execution, or private-code imports.
"""
import json
from pathlib import Path
import sympy as s

a=c=s.Integer(1)
b=s.Rational(5,4);t=s.Rational(3,4);r=s.Rational(1,4);z=c*c
q=1+z*r;h=b-a*a*t;g=(1+z)*b-(1-r+a*a*z*t);j=1+z*r-t
Delta=((1+z)*z*b*b+2*z*(1+z)*b*r+(1+z)*(1+a*a*z)*b
       +z*z*(1-a*a)*r*r+2*z*(1+a*a*z)*r+1+a*a*z)
Nsource=(1+z)*(b+r)*(b+a*a*(1-t))-a*a*(1+a*a*z)*(1-r)*t
source=[s.Integer(1),((b+r)*q*(q+a*a*z*t)-r*Delta)/(2*r*Delta),
        (1+a*a*z)*(1-r)*q*(q+a*a*z*t)/(2*g*Delta),
        q*Nsource/(2*a*a*t*Delta),a*a*q*(b+r)*(q+a*a*z*t)/(2*h*Delta),
        (1+z)*(b+r)*(b-a*a*z*r)*(q+a*a*z*t)/(2*a*a*j*Delta),
        z*z*(b+r)*Nsource/(2*(1-z*h)*Delta)]
assert source==[1,s.Rational(7,66),s.Rational(5,33),s.Rational(5,22),
               s.Rational(10,33),s.Rational(16,33),s.Rational(9,22)]
gaps=[h,g,j,1-z*h,b-a*a*z*r]
assert all(v>0 for v in gaps+source)

u,v=s.symbols('u v',real=True)
W,Xs,Y,Zs=s.symbols('W X Y Z',real=True)
I=s.eye(2);Z=s.diag(1,-1);X=s.Matrix([[0,1],[1,0]])
A0=s.kronecker_product(Z,I)
A1=s.kronecker_product(((1-u*u)*Z+2*u*X)/(1+u*u),I)
B0=s.kronecker_product(I,Z)
B1=s.kronecker_product(I,((1-v*v)*Z+2*v*X)/(1+v*v))
operators=[s.eye(4),A0,A1,B0,B1,A0*B0,A0*B1,A1*B0,A1*B1]
d,n,m,C,A,D,al=source
coefficient=[d+n+m+C+A+D+al,n-A-al,m-C-D,-n+C+A,-m+D-al,
             -d/2-n-A,d/2+al,-d/2-C,-d/2-m-D]
bell=sum((x*y for x,y in zip(coefficient,operators)),s.zeros(4))
psi=s.Matrix([b/a,-c*r,c*h/a,1])
assert psi==s.Matrix([s.Rational(5,4),-s.Rational(1,4),s.Rational(1,2),1])
contact={u:c,v:a*c}
assert bell.subs(contact)*psi==s.zeros(4,1)
assert all(s.factor((psi.T*s.diff(bell,angle).subs(contact)*psi)[0])==0 for angle in(u,v))

# The six-column event-gradient convention in 05_entry.tex EN16.
du=1+u*u;dv=1+v*v
events=[du*dv*Xs*Xs,(v*W-Xs+u*v*Y-u*Zs)**2,dv*(u*W-Y)**2,
        du*dv*Y*Y,(u*W+u*v*Xs-Y-v*Zs)**2,du*(v*Y-Zs)**2]
values={W:psi[0],Xs:psi[1],Y:psi[2],Zs:psi[3],u:c,v:a*c}
stationary=s.Matrix([[s.diff(packet,variable).subs(values) for packet in events]
                     for variable in(W,Xs,Y,Zs,u,v)])
stationary_det=stationary.det()
printed_det=-64*a*c**6*r*g*t*h*j*(1-z*h)*(1+z)**2*(1+a*a*z)**2*Delta
assert stationary_det==printed_det==-297

# Exact seven-root evaluation and the displayed kernel basis.
root_pairs=[(2,0),(4,0),(5,0),(8,0),(3,1),(6,1),(7,1)]
roots=[operators[i]-operators[k] for i,k in root_pairs]
T=s.Matrix.hstack(*[root.subs(contact)*psi for root in roots])
assert T[:,:4].det()==-s.Rational(1,2)
N=s.Matrix([[-1,s.Rational(3,4),s.Rational(9,8)],
            [1,-s.Rational(5,4),-s.Rational(7,8)],
            [3,0,3],[-3,-s.Rational(3,4),-s.Rational(9,8)],
            [1,0,0],[0,1,0],[0,0,1]])
assert T*N==s.zeros(4,3)

# Independent universal reflection algebra: (Alice word, Bob word).
words=[((),()),((0,),()),((1,),()),((),(0,)),((),(1,)),
       ((0,),(0,)),((0,),(1,)),((1,),(0,)),((1,),(1,))]
def reduced(word):
    out=[]
    for letter in word:
        if out and out[-1]==letter:out.pop()
        else:out.append(letter)
    return tuple(out)
def multiply(x,y):
    return reduced(x[0]+y[0]),reduced(x[1]+y[1])
s11,s12,s13,s22,s23,s33=s.symbols('S11 S12 S13 S22 S23 S33')
unknowns=[s11,s12,s13,s22,s23,s33]
S=s.Matrix([[s11,s12,s13],[s12,s22,s23],[s13,s23,s33]])
gram=N*S*N.T
coeff={}
for i in range(7):
    for k in range(7):
        for wi,ci in [(root_pairs[i][0],1),(root_pairs[i][1],-1)]:
            for wk,ck in [(root_pairs[k][0],1),(root_pairs[k][1],-1)]:
                word=multiply(words[wi],words[wk])
                coeff[word]=coeff.get(word,0)+ci*ck*gram[i,k]
forbidden=[((),(0,1)),((0,1),()),((0,1),(0,1)),
           ((0,1),(0,)),((0,1),(1,0)),((0,1),(1,))]
forbidden_names=['B0B1','A0A1','A0A1B0B1','A0A1B0','A0A1B1B0','A0A1B1']
M=s.Matrix([[s.expand(coeff[word]).coeff(variable) for variable in unknowns]
            for word in forbidden])
literal=s.Matrix([[1,s.Rational(7,4),-s.Rational(31,8),0,s.Rational(9,4),-s.Rational(9,8)],
                  [1,-s.Rational(11,4),s.Rational(23,8),-s.Rational(3,2),-3,s.Rational(15,8)],
                  [-9,-s.Rational(9,4),-s.Rational(99,8),0,-s.Rational(9,4),-s.Rational(27,8)],
                  [-3,s.Rational(9,4),-s.Rational(5,8),0,s.Rational(5,4),s.Rational(19,8)],
                  [0,0,0,0,1,0],
                  [3,s.Rational(11,4),s.Rational(33,8),s.Rational(3,2),3,s.Rational(9,8)]])
assert M==literal
assert M.det()==-405
norm=psi.dot(psi)
CHSH=A0*B0-A0*B1+A1*B0+A1*B1
H=s.eye(4)-CHSH/2
Hvalue=s.factor((psi.T*H.subs(contact)*psi)[0]/norm)
pA=s.factor(psi[2]**2/norm)
concurrence=s.factor(2*abs(psi[0]*psi[3]-psi[1]*psi[2])/norm)
assert Hvalue==-s.Rational(17,46) and pA==s.Rational(2,23)
assert concurrence==s.Rational(22,23)
assert 4*A-2*m-2*D==-s.Rational(2,33)
out={'parameters':{'a':'1','c':'1','b':'5/4','t':'3/4','r':'1/4'},
     'active_source_order':['d','n','m','C','A','D','alpha'],
     'source':list(map(str,source)),'strict_gaps_h_g_j_1zh_b_a2zr':list(map(str,gaps)),
     'contact_state':list(map(str,psi)),'kernel_and_two_angular_jets_verified':True,
     'seven_root_action_matrix':[[str(k) for k in T.row(i)] for i in range(4)],
     'first_four_column_minor':'-1/2','displayed_kernel_verified':True,
     'kernel_basis':[[str(k) for k in N.row(i)] for i in range(7)],
     'forbidden_word_order':forbidden_names,'compressed_symmetric_entry_order':list(map(str,unknowns)),
     'forbidden_coefficient_matrix':[[str(k) for k in M.row(i)] for i in range(6)],
     'forbidden_determinant':'-405','unweighted_event_gradient_determinant':'-297',
     'EN16_witness_convention_verified':True,'H_at_contact':str(Hvalue),
     'event_A_at_contact':str(pA),'concurrence':str(concurrence),
     'equal_angle_required_expression_4A_2m_2D':str(4*A-2*m-2*D),
     'generic_scope':'Nonempty relatively Euclidean-open and relative Zariski-open dense subset of the full-dimensional strict native parameter interior; exceptional rational minors may vanish.',
     'extremality_scope':'Canonical extremality extends to full quantum-valid affine polar because H(contact)<0; strict native C2 also has A>0 and p_A>0.',
     'not_checked':'General five-variable stationary determinant expansion; Bell-sector Schur formulas; finite source-symmetry orbit coverage; full C2 Gram admission.'}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
print('All exact source/contact/kernel/universal-word checks passed.')
print('det(T first4)=-1/2; det(forbidden M)=-405; det(event gradients)=-297')
print('H(contact)=',Hvalue,'; concurrence=',concurrence)

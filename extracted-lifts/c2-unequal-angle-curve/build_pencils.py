"""Build the rational curve pencils and check their exact identities.

Python 3 and SymPy; no optimizer, numerical SDP, or parameter sweep.
Only files beside this script are written.
"""
import json
from pathlib import Path
import sympy as s

HERE = Path(__file__).resolve().parent
a, z = s.symbols('a z', real=True)
LOW, HIGH = s.Rational(3, 4), s.Integer(2)
r, sigma = s.Rational(1, 4), s.Rational(1, 8)
e, p = r + sigma, r * sigma
g = a + e - 1
ell = a * (e - p) + p
Q = a * g + ell
Delta = a**3*p + a*a*(e*e-p) + a*p*(2*e-p) + p*p
B = a**3*p + a*a*(e*e+p-e) + a*((3*e-2)*p-p*p) + p*p
source = s.Matrix([
    1, -a*B/(2*(a+1)*Delta),
    (a+sigma)*(1-r)*(a*sigma+r)*ell/(2*(a+1)*g*Delta),
    r*r*(a+sigma)*Q/(2*(1-r)*Delta),
    a*(a+r)*(a+sigma)*ell/(2*(a+1)*Delta),
    (a+r)*(1-sigma)*(a*r+sigma)*ell/(2*(a+1)*g*Delta),
    sigma*sigma*(a+r)*Q/(2*(1-sigma)*Delta)])
k = 2*(a+1)*g*(1-r)*(1-sigma)*Delta
NUM = s.Matrix([
    [-210,-4704,-20286,15120,41664,10752],
    [0,105,-3213,-16968,38304,-5376],
    [126,2457,12285,5544,0,0],
    [-35,56,1435,-9856,3136,14336],
    [0,-105,-2247,-13356,-9408,59136],
    [147,2499,10878,12936,0,0],
    [-15,84,39,-2076,1056,3072]])
C = NUM / s.Integer(262144)
monomials = s.Matrix([a**i for i in range(6)])
assert (k*source-C*monomials).applyfunc(s.cancel) == s.zeros(7, 1)
assert s.cancel(k-21*(a+1)*(8*a-5)*(32*a**3+112*a*a+23*a+1)/131072) == 0
power = s.Poly(s.expand(B.subs(a, LOW+(HIGH-LOW)*z)), z)
bernstein = [s.factor(sum(power.coeff_monomial(z**i)*s.binomial(j,i)/s.binomial(3,i)
                         for i in range(j+1))) for j in range(4)]
assert bernstein == [-s.Rational(497,4096),-s.Rational(1463,6144),
                     -s.Rational(1737,4096),-s.Rational(633,1024)]

# Six source-independent moments and the complete nine-coordinate embedding.
y = s.symbols('y0:6')
H0 = s.Matrix(3, 3, lambda i,j: y[i+j])
H1 = s.Matrix(3, 3, lambda i,j: y[i+j+1])
ML, MR = H1-LOW*H0, HIGH*H0-H1
assert (ML+MR)/(HIGH-LOW) == H0
assert H1.row(0) == H0.row(1) and H1.row(1) == H0.row(2)
source_matrix = s.diag(ML, MR)
embedding = s.zeros(9, 6)
for src, active in enumerate([0,2,3,4,5,6,7]):
    embedding[active,:] = C[src,:]

# Eight visible behavior coordinates, with no contact-dependent features.
ax0,ax1,by0,by1,E00,E01,E10,E11 = s.symbols('a0 a1 b0 b1 E00 E01 E10 E11')
visible = (ax0,ax1,by0,by1,E00,E01,E10,E11)
H = 1-(E00-E01+E10+E11)/2
pn = (1+ax0-by0-E00)/4
pm = (1+ax1-by1-E11)/4
pC = (1-ax1+by0-E10)/4
pA = (1-ax0+by0-E00)/4
pD = (1-ax1+by1-E11)/4
pa = (1-ax0-by1+E01)/4
pairing = s.Matrix([H,4*pn,4*pm,4*pC,4*pA,4*pD,4*pa])
coefficients = (C.T*pairing).applyfunc(s.expand)
poly = s.expand((monomials.T*coefficients)[0])

# Every pair of degree-two interval Grams has exactly this parameterization.
# Five free entries of G, plus H[0,2], leave six unrestricted auxiliaries.
xi = s.symbols('xi1:7')
g00 = (poly.subs(a,HIGH)/(HIGH-LOW)-2*HIGH*xi[0]
       -HIGH**2*(2*xi[1]+xi[2])-2*HIGH**3*xi[3]-HIGH**4*xi[4])
G = s.Matrix([[g00,xi[0],xi[1]],[xi[0],xi[2],xi[3]],[xi[1],xi[3],xi[4]]])
v2 = s.Matrix([1,a,a*a])
numerator = s.Poly(s.expand(poly-(a-LOW)*(v2.T*G*v2)[0]),a)
hpoly, remainder = s.div(numerator,s.Poly(HIGH-a,a))
assert remainder.is_zero and hpoly.degree() <= 4
hc = [hpoly.coeff_monomial(a**i) for i in range(5)]
J = s.Matrix([[hc[0],hc[1]/2,xi[5]],
              [hc[1]/2,hc[2]-2*xi[5],hc[3]/2],
              [xi[5],hc[3]/2,hc[4]]])
behavior_matrix = s.diag(G,J).applyfunc(s.expand)
assert s.expand((a-LOW)*(v2.T*G*v2)[0]+(HIGH-a)*(v2.T*J*v2)[0]-poly) == 0
assert s.expand(s.trace(source_matrix*behavior_matrix)-sum(coefficients[i]*y[i] for i in range(6))) == 0

# The complete parent: source a is independent of testing/state differentiation.
u,v,W,X,Y,Z,t = s.symbols('u v W X Y Z t')
d,n,m,cw,aw,dw,al = source
Du,Dv = 1+u*u,1+v*v
F = (4*n*Du*Dv*X*X+4*m*(v*W-X+u*v*Y-u*Z)**2
     +4*cw*Dv*(u*W-Y)**2+4*aw*Du*Dv*Y*Y
     +(4*dw+2*d)*(u*W+u*v*X-Y-v*Z)**2
     +(4*al+2*d)*Du*(v*Y-Z)**2+2*d*Dv*(X+u*Z)**2
     -2*d*Du*Dv*Z*Z)
contact = {u:t,v:a*t,W:1,X:-t/4,Y:t/4,Z:1}
def at_contact(expression):
    return s.cancel(expression.subs(contact).subs(a,1/(2*t*t)))
assert at_contact(F) == 0
assert all(at_contact(s.diff(F,w)) == 0 for w in (W,X,Y,Z,u,v))

# A strict interval moment point supplies exact dual attainment.
strict = {y[i]:(HIGH**(i+1)-LOW**(i+1))/s.Integer(i+1) for i in range(6)}
strict_minors = [[s.factor(M.subs(strict)[:i,:i].det()) for i in range(1,4)]
                 for M in (ML,MR)]
assert all(v>0 for row in strict_minors for v in row)

def rational(value):
    num,den = s.Rational(value).as_numer_denom()
    return [int(num),int(den)]

def export(matrix,variables):
    records=[]
    rebuilt=s.zeros(matrix.rows)
    for i in range(matrix.rows):
        for j in range(i,matrix.cols):
            p=s.Poly(matrix[i,j],*variables)
            assert p.total_degree()<=1
            values=[p.coeff_monomial(1)]+[p.coeff_monomial(v) for v in variables]
            for col,value in enumerate(values):
                if value:
                    num,den=rational(value)
                    records.append([col,i,j,num,den])
                    term=value*(1 if col==0 else variables[col-1])
                    rebuilt[i,j]+=term
                    if i!=j: rebuilt[j,i]+=term
    assert (rebuilt-matrix).applyfunc(s.expand)==s.zeros(matrix.rows)
    return {'variables':list(map(str,variables)),'matrix_order':6,'blocks':[3,3],
            'record_convention':'[k,i,j,numerator,denominator], zero-based; k=0 constant, k>=1 variable[k-1]. Mirror each off-diagonal record once.',
            'upper_triangle_records':records}

primal=export(source_matrix,y)
primal.update({
    'source_order':['d','lambda','n','m','C','A','D','alpha','beta'],
    'source_embedding':[[rational(v) for v in row] for row in embedding.tolist()],
    'source_equation':'source = source_embedding * (y0,...,y5)^T',
    'active_source_order':['d','n','m','C','A','D','alpha'],
    'active_embedding_numerators':[[int(v) for v in row] for row in NUM.tolist()],
    'active_embedding_denominator':262144,'moment_count':6,
    'parameter_interval':[[3,4],[2,1]],
    'scope':'Exactly the closed cone generated by this fixed unequal-angle curve.'})
dual=export(behavior_matrix,visible+xi)
dual.update({
    'visible_coordinates':list(map(str,visible)),
    'auxiliary_coordinates':list(map(str,xi)),'auxiliary_count':6,
    'polynomial_coefficient_affine':[
        [rational(s.Poly(v,*visible).coeff_monomial(1))]
        +[rational(s.Poly(v,*visible).coeff_monomial(w)) for w in visible]
        for v in coefficients],
    'polynomial_coefficient_affine_order':['constant']+list(map(str,visible)),
    'scope':'Exactly imposes all inequalities from this curve cone on eight-coordinate behaviors; not full Q222.'})
receipt={
    'all_seven_native_curve_rescaling_identities_exact':True,
    'clearing_factor':str(s.factor(k)),
    'clearing_factor_minimum':str(s.factor(k.subs(a,LOW))),
    'B0_Bernstein_coefficients':list(map(str,bernstein)),
    'complete_parent_contact_and_six_first_derivatives_exact':True,
    'source_parameter_held_independent_during_testing_derivatives':True,
    'source_pencil_export_roundtrip_exact':True,
    'behavior_pencil_export_roundtrip_exact':True,
    'coefficient_trace_pairing_exact':True,
    'dual_polynomial_identity_exact':True,
    'strict_moment_leading_principal_minors':[[str(v) for v in row] for row in strict_minors],
    'source_upper_triangle_record_count':len(primal['upper_triangle_records']),
    'behavior_upper_triangle_record_count':len(dual['upper_triangle_records']),
    'source_PSD_order':6,'behavior_PSD_order':6,
    'behavior_auxiliary_count':6,'maximum_representing_atoms':3,
    'no_numerical_optimization_or_parameter_sweep':True}
for filename,data in [('source_pencil.json',primal),('behavior_pencil.json',dual),('exact_checks.json',receipt)]:
    (HERE/filename).write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
print('Exact curve, both pencils, full contact and dual trace identities checked.')
print('PSD orders: source 6, behavior 6; behavior auxiliaries: 6.')
print('Upper-triangle records:',len(primal['upper_triangle_records']),len(dual['upper_triangle_records']))

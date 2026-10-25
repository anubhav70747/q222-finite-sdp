"""Export the exact source and behavior pencils; requires SymPy, no solver."""
import json
from pathlib import Path
import sympy as s

A,m,D,C,alpha=s.symbols('A m D C alpha')
source_vars=(A,m,D,C,alpha)
source=s.diag(s.Matrix([[m,A-m],[A-m,C]]),
              s.Matrix([[D,A-D],[A-D,alpha]]),A-m,A-D,m+D-A)
embedding=s.Matrix([[4,-2,-2,0,0],[0,0,0,0,0],[-1,1,1,0,0],
                   [0,1,0,0,0],[0,0,0,1,0],[1,0,0,0,0],
                   [0,0,1,0,0],[0,0,0,0,1],[0,0,0,0,0]])
a0,a1,b0,b1,e00,e01,e10,e11=s.symbols('a0 a1 b0 b1 E00 E01 E10 E11')
q=(a0,a1,b0,b1,e00,e01,e10,e11)
xi,eta,rho,sigma=s.symbols('xi eta rho sigma')
aux=(xi,eta,rho,sigma)
H=1-(e00-e01+e10+e11)/2
pn=(1+a0-b0-e00)/4
pm=(1+a1-b1-e11)/4
pC=(1-a1+b0-e10)/4
pA=(1-a0+b0-e00)/4
pD=(1-a1+b1-e11)/4
pa=(1-a0-b1+e01)/4
LA=4*(H-pn+pA)
Lm=-2*H+4*(pn+pm)
LD=-2*H+4*(pn+pD)
LC=4*pC
Lalpha=4*pa
dual=s.diag(s.Matrix([[LA+Lm-2*eta-sigma,xi],[xi,LC]]),
            s.Matrix([[LA+LD-2*xi-rho,eta],[eta,Lalpha]]),
            rho,sigma,2*xi+2*eta+rho+sigma-LA)
dual=dual.applyfunc(s.expand)
slack=(4*A-2*m-2*D)*H+4*((m+D-A)*pn+m*pm+C*pC+A*pA+D*pD+alpha*pa)
assert s.expand(s.trace(source*dual)-slack)==0

def export(matrix,variables):
    records=[]
    reconstructed=s.zeros(matrix.rows)
    for i in range(matrix.rows):
        for j in range(i,matrix.cols):
            p=s.Poly(matrix[i,j],*variables)
            assert p.total_degree()<=1
            coeffs=[p.coeff_monomial(1)]+[p.coeff_monomial(v) for v in variables]
            for k,value in enumerate(coeffs):
                if value:
                    num,den=value.as_numer_denom()
                    records.append([k,i,j,int(num),int(den)])
                    term=value*(1 if k==0 else variables[k-1])
                    reconstructed[i,j]+=term
                    if i!=j: reconstructed[j,i]+=term
    assert (reconstructed-matrix).applyfunc(s.expand)==s.zeros(matrix.rows)
    return {'variables':list(map(str,variables)), 'matrix_order':matrix.rows,
            'blocks':[2,2,1,1,1],
            'record_convention':'[k,i,j,numerator,denominator], zero-based indices; k=0 constant, k>=1 indexes variables[k-1]. Mirror each off-diagonal entry once.',
            'upper_triangle_records':records}

primal=export(source,source_vars)
primal['source_order']=['d','lambda','n','m','C','A','D','alpha','beta']
primal['source_embedding']=embedding.tolist()
primal['source_embedding']=[[int(x) for x in row] for row in primal['source_embedding']]
primal['auxiliary_count']=0
behavior=export(dual,q+aux)
behavior['visible_coordinates']=list(map(str,q))
behavior['auxiliary_coordinates']=list(map(str,aux))
behavior['auxiliary_count']=4
behavior['scope']='Exactly enforces all inequalities in the equal-angle C2 source cone. It is an additional quantum-valid constraint, not the full Q222 set.'
for filename,data in [('source_pencil.json',primal),('behavior_pencil.json',behavior)]:
    Path(__file__).with_name(filename).write_text(json.dumps(data,indent=2)+'\n',encoding='utf8',newline='\n')
print('Both exact pencils reconstructed; trace pairing equals the Bell slack.')

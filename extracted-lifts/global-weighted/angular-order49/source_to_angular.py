"""Evaluate the verified quartic source map using exact rational arithmetic.

Usage: python source_to_angular.py d lambda n m C A D alpha beta
Output gives the nine inputs to ANGULAR_CONE_ORDER49.coo.json. This does
not solve the cone problem. Sources outside the nonnegative orthant are
outside the proved Bell-validity equivalence.
"""
from fractions import Fraction
from pathlib import Path
import argparse
import json

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('source',nargs=9,type=Fraction)
args=parser.parse_args()
if any(x<0 for x in args.source):parser.error('all nine canonical source coefficients must be nonnegative')
path=Path(__file__).resolve().parent.parent/'weighted/DETERMINANT_SOURCE_COEFFICIENTS.json'
data=json.loads(path.read_text())
values={}
for row in data['coefficients']:
    value=Fraction(0)
    for exponents,num,den in row['terms']:
        term=Fraction(num,den)
        for x,k in zip(args.source,exponents):term*=x**k
        value+=term
    i,j=row['ij']
    name='h' if (i,j)==(1,1) else f'q_{i-1}_{j-1}'
    values[name]=str(value)
print(json.dumps({'source':dict(zip(data['variables'],map(str,args.source))),
                  'pencil_inputs':values,'map_degree':4},indent=2))

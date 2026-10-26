#!/usr/bin/env python3
"""Check SD2 represents the original unrestricted four-state parent, termwise.
Native coquotients are true source formulas, not independent Bell weights.
"""
from c1_native import *
Lamb=R*T*c3/(2*a)
factor=Delta**2*Lamb
checks={
'm-positive':4*m*factor-4*T*c2*R/M,
'C-positive':4*C*factor-Delta**2*4*R*c3*T/a,
'A-positive':4*A*factor-4*R*T*Delta**2*eta**2*c4,
'D-full':(4*D+2*dB)*factor-Delta**2*4*R*c5*T/J,
'alpha+core-full':(4*al+2*dB)*factor-(4*R*T*Delta**4*k**6*eta**2*cal+2*R*T*Delta**2*eta*cd),
'beta+core-full':(4*be+2*dB)*factor-4*T*cb*R/t,
'signed-core':2*dB*factor-2*R*T*Delta**2*eta*cd,
}
for name,difference in checks.items():
 assert s.cancel(difference)==0,name
 print('PASS',name)
print('PASS: SD2 = Delta^2 Lambda times whole canonical C1 parent, unrestricted W,X,Y,Z and u,v.')

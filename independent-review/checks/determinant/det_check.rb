# Exact rational check of Proof.tex eq. (insight-determinant):
# det P_c = 16 Du^2 Dv^2 u^2 v^2 L_c(u^2,v^2), L_c Laurent with exponents in {-1,0,1}^2, q_e >= 0.
def F(w,u,v,s)
  lam,n,m,c,a,dd,al,be,d = w
  wW,x,y,z = s
  du=1+u*u; dv=1+v*v
  4*lam*du*(wW+v*x)**2 + 4*n*du*dv*x**2 + 4*m*(v*wW-x+u*v*y-u*z)**2 +
  4*c*dv*(u*wW-y)**2 + 4*a*du*dv*y**2 + (4*dd+2*d)*(u*wW+u*v*x-y-v*z)**2 +
  (4*al+2*d)*du*(v*y-z)**2 + (4*be+2*d)*dv*(x+u*z)**2 - 2*d*du*dv*z**2
end
def pmat(w,u,v)
  e=->(i){a=[0,0,0,0];a[i]=1;a}
  p=Array.new(4){Array.new(4,0r)}
  4.times{|i| p[i][i]=F(w,u,v,e.(i))}
  4.times{|i|(i+1...4).each{|j| s=e.(i).zip(e.(j)).map{|a,b|a+b}
    p[i][j]=p[j][i]=(F(w,u,v,s)-p[i][i]-p[j][j])/2}}
  p
end
def det(m)
  m=m.map(&:dup); n=m.size; d=1r
  n.times{|k| piv=(k...n).find{|i| m[i][k]!=0}; return 0r unless piv
    if piv!=k then m[k],m[piv]=m[piv],m[k]; d=-d end
    d*=m[k][k]
    (k+1...n).each{|i| f=m[i][k]/m[k][k]; (k...n).each{|j| m[i][j]-=f*m[k][j]}}}
  d
end
def solve(a,b)
  n=a.size; m=a.each_with_index.map{|r,i| r.dup+[b[i]]}
  n.times{|k| piv=(k...n).find{|i| m[i][k]!=0}; m[k],m[piv]=m[piv],m[k]
    (0...n).each{|i| next if i==k; f=m[i][k]/m[k][k]; (k..n).each{|j| m[i][j]-=f*m[k][j]}}}
  n.times.map{|i| m[i][n]/m[i][i]}
end
EXP=[-1,0,1].product([-1,0,1])
def L(w,u,v); du=1+u*u; dv=1+v*v; det(pmat(w,u,v))/(16*du**2*dv**2*u**2*v**2); end
srand(7)
rr=->{Rational(rand(1..9),rand(1..9))}
names=%w[lam n m C A D alpha beta d]
5.times do |trial|
  w=Array.new(9){rr.()}
  w[0]=0r if trial>=3; w[1]=0r if trial>=3   # C1-like support
  pts=[1,2,3].product([1,2,3]).map{|a,b| [Rational(a+1,a+3),Rational(2*b+1,b+2)]}
  rows=pts.map{|u,v| EXP.map{|e1,e2| (u*u)**e1*(v*v)**e2}}
  coef=solve(rows,pts.map{|u,v| L(w,u,v)})
  # check on fresh points, including odd u,v sign (only even powers should appear)
  ok=6.times.all?{|k| u=Rational(-(k+3),5); v=Rational(k+1,4)
    EXP.zip(coef).sum{|(e1,e2),c| c*(u*u)**e1*(v*v)**e2}==L(w,u,v)}
  qe=EXP.zip(coef).reject{|e,_| e==[0,0]}
  puts "trial #{trial}: Laurent-in-(U,V) on {-1,0,1}^2 exact: #{ok}; min off-centre q_e = #{qe.map{|_,c|c}.min.to_f.round(4)}; all q_e>=0: #{qe.all?{|_,c|c>=0}}; h_c=#{coef[4].to_f.round(3)}"
end

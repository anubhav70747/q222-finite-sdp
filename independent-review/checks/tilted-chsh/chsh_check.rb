# Independent numeric check of the tilted-CHSH pencil (pure Ruby stdlib).
require 'json'
SQ2 = Math.sqrt(2.0)
srand(20261023)

# ---- helpers: symmetric eigenvalues via Jacobi (real symmetric n x n) ----
def jacobi_eigs(a_in)
  n = a_in.size
  a = a_in.map(&:dup)
  200.times do
    off = 0.0
    n.times { |i| n.times { |j| off += a[i][j]**2 if i != j } }
    break if off < 1e-30
    (0...n).each do |p|
      ((p+1)...n).each do |q|
        next if a[p][q].abs < 1e-300
        th = (a[q][q] - a[p][p]) / (2.0 * a[p][q])
        t = (th >= 0 ? 1.0 : -1.0) / (th.abs + Math.sqrt(th*th + 1.0))
        c = 1.0 / Math.sqrt(t*t + 1.0); s = t * c
        n.times do |k|
          akp = a[k][p]; akq = a[k][q]
          a[k][p] = c*akp - s*akq; a[k][q] = s*akp + c*akq
        end
        n.times do |k|
          apk = a[p][k]; aqk = a[q][k]
          a[p][k] = c*apk - s*aqk; a[q][k] = s*apk + c*aqk
        end
      end
    end
  end
  (0...n).map { |i| a[i][i] }.sort
end

# pencil D(t,s) = diag([[sqrt2+t, s/2],[s/2, sqrt2-t]], 1-t, 1+t)  (4x4), as in the repo's DIRECT_ELIMINATED file
def pencil4(t, s)
  [[SQ2 + t, s/2.0, 0, 0], [s/2.0, SQ2 - t, 0, 0], [0, 0, 1 - t, 0], [0, 0, 0, 1 + t]]
end
def pencil_min_eig(t, s)
  jacobi_eigs(pencil4(t, s)).first
end
def in_set(t, s)
  (t.abs <= 1.0 && t*t + s*s/4.0 <= 2.0)
end

puts "== (A) pencil vs set on a grid =="
mism = 0; npts = 0; near = 0
tol = 1e-9
(-1300..1300).step(5) do |ti|
  t = ti / 1000.0
  (-3200..3200).step(8) do |si|
    s = si / 1000.0
    npts += 1
    m = pencil_min_eig(t, s)
    # classify with margin: skip points within tol of the boundary (floating-point ambiguity)
    margin = [1 - t.abs, 2 - t*t - s*s/4.0].min
    if margin.abs < tol then near += 1; next end
    mism += 1 if (m >= 0) != in_set(t, s)
  end
end
puts "grid points=#{npts}, skipped within #{tol} of boundary=#{near}, mismatches=#{mism}"

# boundary probes: points exactly on boundary +/- epsilon along the ellipse and the lines t=+-1
mm = 0; cnt = 0
360.times do |k|
  th = 2*Math::PI*k/360.0
  # ellipse t^2 + s^2/4 = 2 : t = sqrt2 cos th, s = 2 sqrt2 sin th ; only |t|<=1 part is in set
  [1-1e-7, 1+1e-7].each do |f|
    t = SQ2*Math.cos(th)*f; s = 2*SQ2*Math.sin(th)*f
    next if (t.abs - 1).abs < 1e-6
    cnt += 1
    mm += 1 if (pencil_min_eig(t, s) >= 0) != in_set(t, s)
  end
end
[ -1.0, 1.0 ].each do |t|
  # at t=+-1 the set is |s|<=2 (t^2+s^2/4<=2 -> s^2<=4)
  [[0.0,true],[1.999999,true],[2.000001,false],[-1.999999,true],[-2.000001,false]].each do |s,exp|
    cnt += 1
    mm += 1 if (pencil_min_eig(t, s) >= -1e-15) != exp
    mm += 1 if in_set(t, s) != exp && s.abs.round(5) != 2.0
  end
  [1+1e-7].each { |f| cnt += 1; mm += 1 if pencil_min_eig(t*f, 0.0) >= 0 }
end
puts "boundary probes=#{cnt}, mismatches=#{mm}"
# exact-ish check of the determinant identity: det block = 2 - t^2 - s^2/4
maxerr = 0.0
2000.times do
  t = rand*2-1; s = (rand*2-1)*3
  d = (SQ2+t)*(SQ2-t) - s*s/4.0
  maxerr = [maxerr, (d - (2 - t*t - s*s/4.0)).abs].max
end
puts "max |det - (2-t^2-s^2/4)| over 2000 random pts = #{'%.2e' % maxerr}"

# ---- (B) random two-qubit quantum behaviours ----
# 4x4 complex arithmetic by hand using Ruby Complex
Z = [[1,0],[0,-1]]; X = [[0,1],[1,0]]; Y = [[0,Complex(0,-1)],[Complex(0,1),0]]; I2 = [[1,0],[0,1]]
def kron(a, b)
  Array.new(4) { |i| Array.new(4) { |j| a[i/2][j/2] * b[i%2][j%2] } }
end
def expect(psi, m)
  (0...4).sum { |i| (0...4).sum { |j| psi[i].conj * m[i][j] * psi[j] } }.real
end
def obs_xz(theta)
  [[Math.cos(theta), Math.sin(theta)], [Math.sin(theta), -Math.cos(theta)]]  # cos Z + sin X
end
def obs_gen(n)
  Array.new(2) { |i| Array.new(2) { |j| n[0]*X[i][j] + n[1]*Y[i][j] + n[2]*Z[i][j] } }
end
def randn
  Math.sqrt(-2*Math.log(1-rand)) * Math.cos(2*Math::PI*rand)
end
def haar_state
  v = Array.new(4) { Complex(randn, randn) }
  nrm = Math.sqrt(v.sum { |c| c.abs**2 })
  v.map { |c| c / nrm }
end
def unit3
  v = [randn, randn, randn]; n = Math.sqrt(v.sum { |c| c*c }); v.map { |c| c/n }
end
def behaviour(psi, a0, a1, b0, b1)
  ops = { a0: kron(a0, I2), a1: kron(a1, I2), b0: kron(I2, b0), b1: kron(I2, b1) }
  e = { a0: expect(psi, ops[:a0]), a1: expect(psi, ops[:a1]), b0: expect(psi, ops[:b0]), b1: expect(psi, ops[:b1]) }
  e00 = expect(psi, kron(a0, b0)); e01 = expect(psi, kron(a0, b1))
  e10 = expect(psi, kron(a1, b0)); e11 = expect(psi, kron(a1, b1))
  [e[:a0], e[:a1], e[:b0], e[:b1], e00, e01, e10, e11]
end

puts "\n== (B) random quantum behaviours =="
[[:xz, 200_000], [:general, 100_000]].each do |kind, n|
  worst_ell = -Float::INFINITY; worst_pencil = Float::INFINITY; worst_t = 0.0; viol = 0; maxs = 0.0
  n.times do
    psi = haar_state
    if kind == :xz
      a0, a1, b0, b1 = Array.new(4) { obs_xz(rand*2*Math::PI) }
    else
      a0, a1, b0, b1 = Array.new(4) { obs_gen(unit3) }
    end
    q = behaviour(psi, a0, a1, b0, b1)
    t = q[0]; s = q[4] + q[5] + q[6] - q[7]
    ell = t*t + s*s/4.0
    worst_ell = [worst_ell, ell].max
    worst_t = [worst_t, t.abs].max
    mp = pencil_min_eig(t, s)
    worst_pencil = [worst_pencil, mp].min
    viol += 1 if t.abs > 1 + 1e-12 || ell > 2 + 1e-12 || mp < -1e-12
    maxs = [maxs, s.abs].max
  end
  puts "#{kind}: samples=#{n}  max(t^2+s^2/4)=#{'%.12f' % worst_ell} (bound 2)  max|t|=#{'%.12f' % worst_t}  max|s|=#{'%.6f' % maxs} (<=2*sqrt2=#{'%.6f' % (2*SQ2)})  min pencil eig=#{'%.3e' % worst_pencil}  violations(tol 1e-12)=#{viol}"
end

# ---- (C) boundary attainment ----
puts "\n== (C) boundary attainment via the repo's explicit family =="
# |psi_t> = sqrt((1+t)/2)|00> + sqrt((1-t)/2)|11>, A0=Z, A1=X, B_{0,1}=(Z +- r X)/h
maxerr = 0.0; maxres = 0.0; maxrefl = 0.0
(-100..100).each do |k|
  t = k / 100.0
  r = Math.sqrt(1 - t*t); h = Math.sqrt(2 - t*t)
  psi = [Complex(Math.sqrt((1+t)/2)), 0, 0, Complex(Math.sqrt((1-t)/2))]
  b0 = [[1/h, r/h], [r/h, -1/h]]; b1 = [[1/h, -r/h], [-r/h, -1/h]]
  q = behaviour(psi, Z, X, b0, b1)
  exp = [t, 0, t/h, t/h, 1/h, 1/h, (1-t*t)/h, -(1-t*t)/h]
  maxerr = [maxerr, q.zip(exp).map { |u, v| (u - v).abs }.max].max
  s = q[4] + q[5] + q[6] - q[7]
  maxres = [maxres, (s - 2*h).abs, (q[0]**2 + s*s/4.0 - 2).abs].max
  # B observables are reflections (square to identity)
  [b0, b1].each do |b|
    sq = Array.new(2) { |i| Array.new(2) { |j| (0...2).sum { |l| b[i][l]*b[l][j] } } }
    maxrefl = [maxrefl, (sq[0][0]-1).abs, (sq[1][1]-1).abs, sq[0][1].abs].max
  end
end
puts "t grid -1..1 step .01 (201 values): max |behaviour - closed form| = #{'%.2e' % maxerr}; max |s-2h|, |t^2+s^2/4-2| = #{'%.2e' % maxres}; max reflection defect = #{'%.2e' % maxrefl}"

# ---- (D) independent support-function check: max over quantum of alpha*a0 + s  == sqrt(8+2 alpha^2) ----
puts "\n== (D) support function: hill-climb max of alpha*t + s over random 2-qubit models =="
def objective(params, alpha)
  # params: 8 reals for complex state (re,im x4), 4 angles
  v = (0...4).map { |i| Complex(params[2*i], params[2*i+1]) }
  nrm = Math.sqrt(v.sum { |c| c.abs**2 }); v = v.map { |c| c/nrm }
  a0 = obs_xz(params[8]); a1 = obs_xz(params[9]); b0 = obs_xz(params[10]); b1 = obs_xz(params[11])
  q = behaviour(v, a0, a1, b0, b1)
  alpha*q[0] + q[4] + q[5] + q[6] - q[7]
end
[-3.0, -2.0, -1.0, 0.0, 0.5, 1.0, 1.7, 2.0, 2.5].each do |alpha|
  best = -Float::INFINITY
  6.times do
    p = Array.new(8) { randn } + Array.new(4) { rand*2*Math::PI }
    f = objective(p, alpha); step = 0.5
    4000.times do
      cand = p.map { |x| x + step*randn*0.5 }
      fc = objective(cand, alpha)
      if fc > f then p, f = cand, fc else step *= 0.999 end
      step = [step, 1e-7].max
    end
    best = [best, f].max
  end
  bound = Math.sqrt(8 + 2*alpha*alpha)
  # alpha=|alpha|>2 case in the repo: support is 2+|alpha| only for the OUTER tilted description; true quantum value is sqrt(8+2 alpha^2) for |alpha|<=2 ... report both
  puts "alpha=%5.2f  hill-climb best=%.9f   sqrt(8+2a^2)=%.9f   gap=%.2e   (2+|a|=%.6f)" % [alpha, best, bound, bound-best, 2+alpha.abs]
end

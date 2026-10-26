"""run_points: run the exact primal admission certificate (exact_primal.py) on a fixed list of sources and a negative control,
print a results table and write admission_points.json next to this script.

usage:  python run_points.py [--quiet] [--skip-control] [--out FILE]
Every listed point must give exact PD True with 0 residual rows; the out-of-chart control must NOT certify.  Exit status 0 iff so."""
import argparse, json, os, platform, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from exact_primal import exact_primal

POINTS = [('C1', (0.111, 0.0967, 0.0543, 0.114, 1.0)),
          ('C1', (0.3645, 0.01901, 0.01046, 0.8042, 1.086)),
          ('C1', (0.0945, 0.146, 0.0544, 0.335, 0.999)),
          ('C2', (0.02828, 0.0009625, 0.001251, 0.47, 0.9988)),
          ('C2', (1.0, 0.5, 0.65, 0.4, 0.3))]
CONTROL = ('C2', (1.0, 0.5, 0.65, 0.4, 0.6))        # n < 0: outside the native chart, so it must NOT certify


def run_one(fam, p, verbose, allow_outside=False):
    t0 = time.time()
    try: r = exact_primal(fam, list(p), verbose=verbose, allow_outside=allow_outside)
    except Exception as ex: r = dict(error='%s: %s' % (type(ex).__name__, ex))
    rec = dict(family=fam, point=list(p), allow_outside=allow_outside, rational_point=r.get('th'), t=r.get('t'), trel=r.get('trel'), solver_status=r.get('status'),
               exact_residual_rows=r.get('exact_residual_rows'), exact_pd_blocks=r.get('exact_pd_blocks'), exact_pd=bool(r.get('exact_pd', False)),
               time_s=round(time.time() - t0, 1))
    if 'error' in r: rec['error'] = r['error']
    return rec


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--quiet', action='store_true', help='suppress the per-point trace of exact_primal')
    ap.add_argument('--skip-control', action='store_true', help='do not run the negative control')
    ap.add_argument('--out', default=os.path.join(HERE, 'admission_points.json'))
    a = ap.parse_args(); T0 = time.time(); recs = []
    for fam, p in POINTS:
        print('\n=== %s %s' % (fam, ' '.join(map(str, p))), flush=True); rec = run_one(fam, p, not a.quiet); rec['expected_exact_pd'] = True; recs.append(rec)
    if not a.skip_control:
        fam, p = CONTROL
        print('\n=== negative control %s %s  (default: refused as outside the native chart)' % (fam, ' '.join(map(str, p))), flush=True)
        rec = run_one(fam, p, not a.quiet); rec['expected_exact_pd'] = False; recs.append(rec)
        print('\n=== negative control %s %s  with --allow-outside semantics' % (fam, ' '.join(map(str, p))), flush=True)
        rec = run_one(fam, p, not a.quiet, allow_outside=True); rec['expected_exact_pd'] = False; recs.append(rec)
    for r in recs: r['as_expected'] = (r['exact_pd'] == r['expected_exact_pd']) and (not r['expected_exact_pd'] or r['exact_residual_rows'] == 0)
    total = time.time() - T0; ok = all(r['as_expected'] for r in recs)
    import numpy, scipy, clarabel, flint, sympy
    versions = dict(python=platform.python_version(), numpy=numpy.__version__, scipy=scipy.__version__, sympy=sympy.__version__, python_flint=flint.__version__,
               clarabel=getattr(clarabel, '__version__', None))
    json.dump(dict(versions=versions, total_time_s=round(total, 1), all_as_expected=ok, results=recs), open(a.out, 'w'), indent=1); open(a.out, 'a').write('\n')
    W = '%-50s %-11s %-10s %-13s %-8s %-26s %-8s %-3s %s'
    print('\n' + W % ('point', 't', 'trel', 'solver', 'res rows', 'exact PD per block', 'exact PD', 'ok', 'time'))
    for r in recs:
        name = '%s %s%s' % (r['family'], ' '.join('%g' % x for x in r['point']), ' [control, --allow-outside]' if r['allow_outside'] else (' [control]' if not r['expected_exact_pd'] else ''))
        f = lambda v, fmt: fmt % v if v is not None else '-'
        print(W % (name, f(r['t'], '%.4e'), f(r['trel'], '%.3e'), r['solver_status'] or '-', f(r['exact_residual_rows'], '%d'),
                   r['exact_pd_blocks'] if r['exact_pd_blocks'] is not None else r.get('error', '-')[:26], r['exact_pd'], 'yes' if r['as_expected'] else 'NO', '%.0fs' % r['time_s']))
    print('\nall as expected: %s   total %.0fs   wrote %s' % (ok, total, a.out)); return 0 if ok else 1


if __name__ == '__main__': sys.exit(main())

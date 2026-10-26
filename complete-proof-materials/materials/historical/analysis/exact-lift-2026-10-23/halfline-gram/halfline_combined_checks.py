"""Check literal selector/B matrices in the staged TeX against frozen data.

Also records exact witness replay and compiler results already produced by
their separate standard-library checkers. Writes one local check report.
"""
from pathlib import Path
from fractions import Fraction as Q
import hashlib
import json
import re

HERE = Path(__file__).resolve().parent
tex = (HERE/'halfline_combined.tex').read_text()
path = HERE/'SIMPLE_SOURCE_HALFLINE_EXACT_WITNESS.json'
raw = path.read_bytes()
witness = json.loads(raw)

def displayed_matrix(prefix):
    match = re.search(re.escape(prefix)+r'\s*\\begin\{pmatrix\}(.*?)\\end\{pmatrix\}', tex, re.S)
    assert match, prefix
    return [[entry.strip() for entry in row.split('&')]
            for row in match.group(1).strip().split(r'\\')]

for name, prefix in [('BG', r'B_G&=\frac1{256}'), ('BH', r'B_H&=\frac1{18496}')]:
    displayed = [[int(entry) for entry in row] for row in displayed_matrix(prefix)]
    assert displayed == [[int(entry) for entry in row] for row in witness[name+'_compact']['integer_numerator']]
    assert witness[name+'_compact']['denominator'] == ('256' if name == 'BG' else '18496')

j = displayed_matrix(r'J_+(v)=\frac1{64}')
coefficients = [[[0, 0], [0, 0]] for _ in range(3)]
for r in range(2):
    for s in range(2):
        expression = re.sub(r'\s+', '', j[r][s])
        tokens = re.findall(r'[+-]?\d+(?:v(?:\^2)?)?', expression)
        assert ''.join(tokens) == expression
        for token in tokens:
            match = re.fullmatch(r'([+-]?\d+)(v(?:\^2)?)?', token)
            degree = {None: 0, 'v': 1, 'v^2': 2}[match.group(2)]
            coefficients[degree][r][s] += int(match.group(1))
assert coefficients == witness['selector_integer_coefficients']
assert list(map(Q, witness['source'])) == list(map(Q, ['359','0','0','75/16','60','416','41/2','831/2','151/2']))
assert r'\frac{75}{16},60,416' in tex
assert r'\frac{41}{2},\frac{831}{2},\frac{151}{2}' in tex

replay = json.loads((HERE/'simple_source_exact_witness_checks.json').read_text())
assert replay['witness_sha256'] == hashlib.sha256(raw).hexdigest()
assert replay['all_four_blocks_psd']
compiler = json.loads((HERE/'halfline_gram_checks.json').read_text())
assert compiler['auxiliary_count'] == 44 and compiler['canonical_order'] == 33
assert compiler['feasible_auxiliary_supnorm_bound_constant'] == 16
dual = json.loads((HERE/'strict_trace_certificate_checks.json').read_text())
assert dual['input_sha256'] == hashlib.sha256((HERE/'HALFLINE_QUADRATIC_SELECTOR_ORDER24.coo.json').read_bytes()).hexdigest()
assert list(map(Q, dual['source_trace_coefficients'])) == [Q(16)]*9
assert list(map(Q, dual['auxiliary_trace_coefficients'])) == [Q(0)]*44
assert Q(dual['W_lower_eigenvalue_exact']) == Q(1, 2)
log = (HERE/'halfline_combined.log').read_text()
assert not any(token in log for token in ['Overfull', 'undefined references', 'undefined on input'])
match = re.search(r'Output written on halfline_combined.pdf \((\d+) pages,', log)
assert match and int(match.group(1)) == 4
out = {
    'staged_pdf_pages': 4, 'literal_J_BG_BH_match_frozen_witness': True,
    'frozen_witness_sha256': hashlib.sha256(raw).hexdigest(),
    'exact_witness_replay': replay['all_four_blocks_psd'],
    'gram_ranks': {name: row['rank'] for name, row in replay['exact_blocks'].items()},
    'compiler_auxiliaries': 44, 'compiler_orders': [24, 33],
    'uniform_auxiliary_bound_constant': 16, 'latex_references_resolved': True,
    'independent_strict_dual_trace_certificate': True,
    'dual_weight_lower_eigenvalue': '1/2',
    'overfull_boxes': False,
    'scope': 'Private staging of exact quadratic half-line selector cone and one rational rescue; universal adequacy remains open.',
}
(HERE/'halfline_combined_checks.json').write_text(json.dumps(out, indent=2)+'\n')
print(json.dumps(out, indent=2))

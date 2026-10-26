"""Build the human audit's explicit component/leaf records from pinned originals."""
import csv
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent.parent
INDEX=json.loads((ROOT/'recovery/full_r21/records/COMPONENT_INDEX.json').read_text(encoding='utf-8'))
INDEX={x['component_id']:x for x in INDEX['components']}
rows=[]

def add(ids, role, hypothesis, conclusion, dependencies, mechanism, checks='', status='inspected complete operative derivation'):
    for component in ids.split():
        item=INDEX[component]
        source=ROOT/'recovery/full_r21/components/source'/f'{component}.tex'
        rows.append(dict(component_id=component,title=item['title'],original_source_path=item['source'],
            inspected_lines=f'1-{len(source.read_text(encoding="utf-8").splitlines())}',
            role=role,hypotheses=hypothesis,conclusion=conclusion,operative_dependencies=dependencies,
            verified_mechanism=mechanism,fresh_exact_checks=checks,status=status))

add('pe1088004af65 pe351eb06b28b pec6c617d3a3d p49f367f1b5e8','assembly and exhaustive source routing',
    'Original physical finite-contact source over arbitrary real closed Kc; C2 lambda=beta=0; d>0; strict EN gaps; n>=0; independent Kt. CT: complementary C/alpha strong or weak; HJ: directed strong-C weak-alpha actual II with weak total scale.',
    'CT routes every actual-angle type; HJ routes all directed weak-total II sources into WA/WE/SL/WC or R1-R10; LC conditional assembly then universal tensor SOS.',
    'EN/SE/AU/A6 and C1 shared; strong C2; two-weak radial; four regular exits; ten HJ rows',
    'Positive source identities; exact C2 reflection including homogeneous factor; ordered I/U/O partition; strong-C excludes simultaneous r,g,h,j infinitesimal; bare-source BSR hypotheses matched to actual local bodies.',
    '', 'complete routing inspected; global inputs shared with lead/C1 reviewer')

add('pfb3a83f3849a','strong II exact source and full parent reduction',
    'C and alpha have bounded reciprocals; actual u_c,v_c infinitesimal; all six event scales arbitrary.',
    'Full parent in contact graph normals and angular contact ideal with opening-safe bounds.',
    'Original C2 entry; LII local and remote',
    'DI1-DI8 invertible polynomial state graph; exact normal block and both derivatives; bounded products m/z^2,A/y,CL,Dyj,alpha(1-H).')
add('pe16499ccc4c6','strong II local',
    'Strong II preceding source; independently small testing u,v; arbitrary state.',
    'Entire parent tensor SOS at local tests.',
    'DI source; strict bounded Gram',
    'LII4 harmonic completion; LII8-LII15 finite monic full scalar/normal division; exact LII20 zero supplies T pivot; xi/R Schur 2/[h(1+h)^2]>0. No testing contact identified with source contact.')
add('p590112f102bd','strong II remote',
    'Same strong II source; all tests outside the local contact neighborhood.',
    'Every remote finite, reciprocal and literal endpoint test.',
    'DI; strict bounded Gram; A6 shared for noncommuting residues',
    'RR3 complete homogeneous parent; RR5-RR12 normalized harmonic frame and nonempty real allocation interval; b-inf exact three-square RR13; no invented pivot along endpoint null state.')

add('pc3daa35e9a9b pa142a027a282 p76545700d90b p6b22e944bd7c','strong UI full-source/local/all-tests chain',
    'Actual u_c a bounded unit, v_c infinitesimal; C/alpha strong; arbitrary positive complementary openings and all event upper scales.',
    'Entire parent tensor SOS at every independent test/state.',
    'EN; strict Gram; A6 shared',
    'UR canceled products and entire centered packets; reciprocal-contact real matrix identity; UI1 exact zero redistributes soft pivots; UI2 full scalar and normal monic division; positive matched Hessian; UA opening-safe real interior approximation and finite real caps; all endpoint frames.')

add('p54b8f20ad2dd pcee1372511ed p888d9f79ba3b p4363feef60dd','strong HI source/local/remote chain',
    'Actual first angle unbounded, second infinitesimal; both complementary events strong; arbitrary n>=0 including unbounded n and n=0.',
    'Full parent tensor SOS for all tests.',
    'WE3 source exclusion; generic IC/T3 circuit lift shared; strict Gram; A6',
    'Source r/q bounded from positive physical numerator; entire four toric packets retained; N0M arrays and half-angle conductors; full-HI retains actual centered n packet and changes bounded toric coefficients by nu=2nR; pivots lambda(r)=2(1+r)^2/(3+r)>0 and 4/(3+r)>0; remote real caps on complete kernel forms.')

add('pf5726218aec8 pf3e0ca8bf054 p7405dbf89341','strong IH source/local/remote',
    'Actual first angle infinitesimal, second unbounded; C/alpha strong; every original event opening retained.',
    'All independent finite and projective tests/states.',
    'IC exact parent identities; endpoint WIE2 zero; strict Gram; A6',
    'Strong source B unit, r*y bounded; cap C in exact parent without altering source graph; SL4-SL5 positive determinants; complete full normals; all remote angle types including opening kernels and double infinity.')
add('p03eb42346077 p93b05ea6ad1e pc00b84929b12 p18dbedfff317 p0fb32cd2f1bc pf95b1c5fc32f','weak IH all-tests chain',
    'Actual first angle infinitesimal, second unbounded; C weak (alpha strong); B unbounded; rho=B/y bounded in [0,1] in residue; n>=0.',
    'All tests through rho=0/interior/1 contact modules, open-axis and remote modules.',
    'Original EN source; strict Gram; A6',
    'Exact homogeneous IC3 parent; scalar IC7-8 with exact zero; complete IC10-12 normal quotients; rho endpoints retain full D row and normalized X/Y; WIE2/WIE4 exact repairs; remote extra sqrt(m) root; OA7-9 explicit nonempty allocation interval. Rational-curve example not a general theorem premise.')

add('p534f57e1da6f p8b53b68d8af6 p9f1235511111 p53752031d16c p163e0ebf59f6 pc7dacaa3258e','two-weak complementary events',
    'Original physical C2 parent, d=1 and C,alpha infinitesimal; actual contact arbitrary.',
    'Chart-free source entry then complete all-tests tensor SOS.',
    'Generic monic contact algebra; strict Gram; A6',
    'Matching amplitudes and four-gradient Jacobian determine source signs; ordinary coalescence uses all eight lower-bounded matching strengths; small-angle radial source has all gaps infinitesimal. Four complete current packets; exact fourth-current zero; weighted N=sum(gaps) scalar/normal quotients; TN=sqrt(N)T supplies unit radial pivot; explicit remote homogeneous identities. Unit-N optional chart not needed for actual two-weak radial entry.')
add('pf8e9d3c47c72','matching-strong positive-chart predecessor',
    'Physical positive rational C2 chart; C/alpha infinitesimal; matching n,m,A,D have bounded reciprocals.',
    'All tests after source entry into ordinary or small-angle radial chain.',
    'Matching budget; original source denominator; weak-event entry; radial theorems',
    'Budget makes matching masses unbounded; normalized contact matching amplitudes force equal ordinary or paired endpoint angles; exact reflection; a^2 unbounded would force bounded b and infinitesimal t hence infinitesimal D contradiction; bounded-a weak-event entry applies.',
    '', 'complete proof inspected; redundant predecessor superseded by chart-free two-weak entry in final CT route')

add('p7294aae40ee2 p81893ad08bed pd7bd7f0d7fb6 p8f912583feaa pec59340d644f','single-angle finite-b source openings',
    'Actual first angle infinitesimal, second bounded unit; finite native b; r/g/b-h openings arbitrary, including n=0.',
    'All tests for these finite-b source cases.',
    'Generic contact monic division; strict Gram; A6',
    'Source lower b bound from n>=0; canceled C/A products; normalized entire Y row; triangular centered m row at g small; joint n/m entire excess when both open; exact balanced Laurent reciprocal-pair plus three-cycle identity; strict matched conductors and all outer frames.')
add('pb60046a37191 pbb29c73cc124 pab6222b899ef','single-angle unbounded-b small-B',
    'Actual first angle infinitesimal, second unit; native b unbounded; B=x*b infinitesimal; sigma=1/b infinitesimal; no rate relation B/sigma.',
    'Whole parent all-tests tensor SOS.',
    'HL homogeneous scalar lemma; complete normal inversion; strict Gram; A6',
    'Exact canceled source scales and HL scalar support; zero HL11 supplies T pivot at zero sigma residue; UB source normal arrays preserve n=0; full X-normal coupling retained; exact UB17 zero; harmonic double-zero state frame; positive real outer caps at source 1/B scale.')
add('p0c7a15f800d9 p6c62a53b1454 pf7943b34a7dc','strong complementary single-angle pitchfork',
    'Actual first angle infinitesimal, second unit; d positive unit; C,alpha bounded below with arbitrary upper scales; all n,m,A,D>=0.',
    'Entire parent all independent tests/states.',
    'Original Laurent stationary dictionary; strict Gram; A6',
    'Physical test forces A,D and m/k^4 lower bounds; complete graph and both first-normal derivatives; exact finite Taylor scalar division; real hard Schur quartic has double root at nonzero contact and positive leading coefficient; SC9-SC19 covers unbounded C/alpha while retaining each full large packet; SC21-SC23 strict soft pivot; explicit endpoint coordinate Grams.')

add('pe56fb394b522 pb7189b9c1f6e p80253fe9bc1b','angle exclusions and ordinary source orbit',
    'Original C2 physical contact and n>=0; directed weak-C or ordinary actual angles as stated.',
    'Excluded directed weak-C II/UI/IU/HI/HU cells; ordinary-ordinary either two-weak or five strong packets; five consecutive lower-bounded packets all-tests.',
    'Original source; two-weak theorem; A6; strict Gram',
    'Positive numerator implications; WE2-WE6 prove r/q bounded then C reciprocal bounded; five-event finite caps with open-axis kernels, three exact five-coordinate corner quadratics and all higher monomials allocated. Historical example at final angle-reductions paragraph not a coverage premise.')

add('pbaab204a8ab1 p9e725a065583 p6db37c4233da p4b75bfe0436c p261f71af685d pf77de02d80c4','HJ R1 complete source row',
    'B,z,y,t infinitesimal; original physical C2 source, strict gaps and n>=0; native ratio a arbitrary; source j*T arbitrary.',
    'Whole original R1 parent all tests/state signs/endpoints.',
    'Original EN source; strict Gram',
    'R1N exact source map and bounded canceled costs remove former jT restriction; RSR exact scalar identity and three exhaustive signed-cost cases; RFN completes correlated n packet BEFORE scalar spending, full two normals RFN21-24 and exact radial helper Gram; SJ full core Schur/sign reduction; R1N13-24 outer current identities and four finite cap cases.',
    'RFN18c/RFN19 exact polynomial identities and positive principal minors; source-independent current budgets reconstructed afresh.')

add('p913a09bcdcda p393d3faa269a p09b2107cc05e p667905169990 p8f8780350b13 pb88ccd68a4bf pef021eabacde','HJ R2 complete source row',
    'B,z,y infinitesimal, t unbounded; original physical C2 source; a arbitrary and n>=0.',
    'Whole R2 parent all tests and arbitrary state ratio.',
    'R1 generic radial/circuit budgets; SC/QI algebraic conductors; strict Gram',
    'RER exact physical inequality bounds a*y and a(a+1)y before extraction; RL harmonic shift and complete normal budget for a bounded/p small; TU/UL relative current split retains exact extra n normal packet; finite matched circuits with unbounded current retained as rows; R2O exact full Bell reflection then two-current/six-root real cap joins; SJ state signs/endpoints.',
    'SH1 parameterized opposite budget independently expanded; R1 generic budgets independently expanded.')
add('p95350f69a353 pc376284fc7ae','shared SC/QI conductor bodies actually inspected in prior focused C2 audit',
    'Their native source charts for stated theorems; in R2 only explicit polynomial/conductor identities with TU/UL coefficient hypotheses are imported.',
    'Finite exact contact scalar and both normal lift on prescribed root lists.',
    'Monic divisions; strict Gram',
    'SC16-SC20 toric and half-root conductors, three gradient differences with unique relation; QI7-QI8 extra full-normal roots and bounded source coefficient arrays. No C1 physical theorem transported to negative beta.',
    '', 'prior focused complete mechanism inspection; current import hypotheses rechecked')

add('p1c9902e9a6f2 p78dfeedee062 p1de8b5ffab69','HJ R3 and R4 complete rows',
    'R3 t,h units l infinitesimal, or R4 t,h,l units j infinitesimal; e=1-r arbitrary positive; all physical source signs retained.',
    'All tests through e-unit and e-infinitesimal alternatives.',
    'WA/WE full Gram; R5 centered m; RD D gauge; RFN/RSR compact source; strict Gram',
    'R3E/R4E normalize nr/mg/Cl/Dj before estimates and use WE exact a quotient; D row jointly normalized; full local strict 2x2 blocks; centered m at g small; EN small-e exact full normal inverse and real-cap outer frames; complete homogeneous R4E14-16 endpoints.',
    'All15 RD4 entries and all15 centered-m matrix entries independently checked.')

add('p97700cec8dc8 p9ec3d6bbd717 p0d4eedd6c7a4 p296e824cb97a p99372d145f17 p34609a6519a6','regular exits WA WE SL',
    't,g,j units; h infinitesimal; WA r,e,a units; WE r,e units a infinitesimal; SL r infinitesimal and a bounded; actual source and gaps.',
    'Full parent all tests on these regular exits.',
    'WA five-root source-contact identities and strict Gram',
    'WA5-WA8 retain scalar, both full normal rows and hard block; WE monic (b2-alpha*)/a identity regular jointly through h/a arbitrary; SL separate n/A normalizations with no n/A ratio; complete finite residue determinants and reciprocal root ideals. p9ec and p977 differ only in the historical extension reference, not proof content.')
add('p0c68f3aecc59 p7c5e0b890910 pbb649976abc6','regular exit WC',
    't,g,j,a units; h,e infinitesimal; original physical source includes n=0.',
    'Entire WC parent all tests.',
    'IC/T3 restricted contact lift shared; strict Gram; WA source arrays',
    'BS exact positive physical numerator bounds B/e,n/e; WC retains entire qs=nr^2+(m+1/2)e^2; exact monic B divisions of scalar and complete normal arrays; bounded test-ratio conductors and reciprocal exact full-state frame with n/m bounded; GR2-7 six-/five-root remote parents and literal endpoints.')

add('p8da9ea6ed9e5','HJ R5',
    't,h,l,j units and r infinitesimal; a bounded by g>0, residue rho in [0,1].',
    'Whole R5 parent all tests.',
    'WA/WE; strict Gram',
    'Three exhaustive rho cases: interior bounded m, rho=0 scaled R2, rho=1 entire centered m excess. Positive eta/R1 determinant; complete remote determinants, rho=0 extra aW root, single/double reciprocal forms.',
    'Centered-m whole-matrix identity independently checked entrywise.')
add('pd5fa361b4623 pc82cd677f4b9 pb60eba754641 p75b395956985','HJ R6 complete row',
    't,h,l,j,r units, g infinitesimal. B,z,y infinitesimal. a bounded or unbounded; in latter b=a^2B forced bounded and split unit/infinitesimal. n=0 retained where possible.',
    'Every R6 source all tests.',
    'R5 centered m; generic CC/SE conductor shared; BE selected state identities; strict Gram; R6 remote kernels',
    'Bounded-a strict centered-excess block. RB unit-b retains full m and centered n packets, bounded entire inverse/scalar/normals. RI infinitesimal-b proves exact F0/b,AN/b cancellations, restricted five-term raw support, explicit scaled conductor RI14-15, coefficientwise filtrations and positive final 2x2 Schur. Original n used on unit-state chart. RT thin proof consistent but redundant after RB+RI split; no source theorem used at negative beta.',
    'Actual source jets RI18; each raw angular coefficient low weighted part RB/RI; full M normal low coefficients; F0/AN exact face cancellations and leading pivots independently checked.')
add('p9a223c96c696 p405d96527133','HJ R7 complete row',
    't unit h,g infinitesimal; r,j not both infinitesimal; a bounded or unbounded; all original signs/gaps.',
    'Whole R7 row all tests.',
    'R5 centered m; RD exact D gauge; RI scaled conductor; CC shared; BE selected state frames; strict Gram',
    'Bounded-a branch jointly removes entire m/D packets and normalizes divergent n/A; all real common-kernel remote forms strict. Unbounded-a RU5-RU16 retains four entire packets, proves bounded raw scalar and both normals with no h/j pole; RU20 coefficientwise two gradings and exact bounded conductor; full residue determinant (t+2)/(4t)>0 even t=1; all capped remote kernels/homogeneous frames.',
    'RU22 source jets; each raw angular coefficient unit-/small-b low filtration and M normal coefficients independently reconstructed.')
add('p482fcad1a2e2 pe07041e033da','HJ R8 complete row',
    't,g units, h,j infinitesimal; r unit or infinitesimal; bounded a/e openings or unbounded a; all source signs including n=0.',
    'Whole R8 parent all tests.',
    'RD/WA/WE exact Gram; R10 coupled identity/thin stress; OP4 actual n=0 C1 import; SEO/BE joins; CA/SE conductor shared; strict Gram',
    'Small-r jointly normalized sqrt(n),sqrt(A),sqrt(D) five-root Gram and exact D syzygy. Unit-r JD4 removes whole divergent D block; bounded-a positive Schurs; unbounded thick/thin/n=0 split; JD15 paired-D identity cancels j inverse in stress; JD17 raw support/normal filtration; all finite remote and projective frames.',
    'All15 RD4 entries freshly checked against independently reconstructed WA matrix.')
add('pd709027c7805 pb18dc55f4f2f pe2fc3cbe3819','HJ R9 complete row and remote extension',
    't,g,j units; h,e,a infinitesimal; original n>=0 and all strict source gaps.',
    'Whole R9 parent all tests; remote frontier comments on then-open local problem superseded by RC body.',
    'Original source arrays; GR complete remote identities; strict Gram',
    'Canceled source circuit coefficients lambda~1/2 and cF>0; RC coupled full rank packet retains both normals and unbounded testing ratio inside positive roots; final E/RD determinant positive; homogeneous roots remove ratio inverse; frontier WF5/WF6 establishes exact GR remote bounds without a-unit assumption.',
    'All six exact candidate-source stationary equations independently rebuilt and verified modulo k^2=aB,l^2=tB.', 'RC prior complete mechanism inspection; source and remote extension newly checked')
add('p26da50239f47','R9 bounded-testing-ratio partial predecessor',
    'Physical R9 source; local u,v nonzero; bounded pure-testing v/u.',
    'Full local parent on bounded-ratio tests only.',
    'WC complete scalar/normal arrays; R9 source regularity; strict Gram',
    'R9B1 turns the negative circuit into an infinitesimal off-diagonal; R9B4 allocates all complete normal monomials directly; no contact conductor needed. Later RC pb18dc55f4f2f covers every testing ratio and does not import R9B.',
    '', 'full body inspected; redundant partial predecessor, no importing coefficient obligation in final RC route')
add('p526dc49d5fcb p9ee4ba7326a7 pef37475d831e pf1aac1dff7c7','HJ R10 complete row',
    't,g,j units; h,e infinitesimal; a unbounded. b=a^2B and w=ae bounded, b/w bounded by exact n>=0 identity; rho=n/e^2 any nonnegative type.',
    'Entire R10 source all tests/states.',
    'RC coupled identities; OP4 n=0 physical C1 import; CA/SE conductors shared; SEO critical wall; GR remote; strict Gram',
    'Full weighted parent R10.9 and strict thick rho Schur; literal n=0 exact C1 source; positive thin layer TL3 full coefficient comparison with exact delta m=0, new centered n row TL4, entire stress TL6 and both normal/raw filters TL8-10; formal beta used only in positive paired packets. UW polynomial numerator divisible by s*w gives deficit/m=s bounded and full unit-W correction. Maximal state partition supplies exactly the contact chart.',
    'SEO full determinant, zero section and outer identities independently reconstructed.')
add('pc024f5300388','literal zero overlap and exact formal comparison',
    'Original weak C2 source; n>=0 including equality.',
    'OP4 same actual contact; physical C1 source only at n=0; positive n yields beta1<0 and is used algebraically only.',
    'Original EN; C1 source uniqueness shared',
    'Exact OP2 positive source partition and OP5 beta=-sn/2. All C2 importing maps explicitly retain native eta=t, zeta=1/a, r=h and e_scale^2=t*a^2B; correct C1 theorem selected by scale.')

add('p05a86be0046c p69995d829bb2','shared critical-wall and selected state/projective identities',
    'In C2 imports only: d positive unit; C,A,D have positive lower bounds. SEO12 case C,d bounded units, A,D may be arbitrarily large. BE4 requires m,A lower bounds and Z/W infinitesimal. BE infinity selected forms require only displayed C,A or C,A,D lower bounds.',
    'Complete selected parent strict Grams at all imported remote/unit-state/endpoints.',
    'Strict Gram; original event packets; pure-testing state normalization',
    'SEO11 exact full parent and opening-safe normalized A/D roots; SEO13 strict uV/vV scalar masses; complete hard-V rows degree>=2; SEO17/18 outer identities and SEO20 determinant; BE4 entire (X,Y,uW,vW) root ideal; BE2/infinity0 entire homogeneous root ideals, not state restrictions. Actual C2 import source bounds checked in each leaf.',
    'SEO13/17/18/20 exact polynomial identities checked independently; critical hard restriction visibly strict.')
add('p25af2c58da7e','strict bounded tensor Gram lemma',
    'Finite symmetric Gram over Bc tensor Bt has positive definite real residue; arbitrary tensor row entries.',
    'Finite SOS in Kc tensor_R Kt; no mixed inverse or radical.',
    'Real spectral/Cholesky decomposition of pure-field matrices',
    'Finite positive-real budget splits each bounded pure-field product perturbation into separated nonnegative factors; whiten only a pure-field matrix. No order on tensor ring and no algebraic disjointness hypothesis.',
    '', 'prior focused complete proof inspection')

with (HERE/'component_coverage.csv').open('w',encoding='utf-8',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=list(rows[0]))
    writer.writeheader(); writer.writerows(rows)
(HERE/'component_coverage.json').write_text(json.dumps(rows,indent=2,ensure_ascii=False),encoding='utf-8')
(HERE/'inspected_component_ids.json').write_text(json.dumps([r['component_id'] for r in rows],indent=2),encoding='utf-8')
print(json.dumps({'components':len(rows),'unique':len({r['component_id'] for r in rows}),'all_have_source':True}))

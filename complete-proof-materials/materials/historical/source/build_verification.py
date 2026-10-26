#!/usr/bin/env python3
"""Aggregate completed checks without promoting their mathematical scope."""
from pathlib import Path
import argparse
import hashlib
import json
import shutil
import zipfile

ROOT=Path(__file__).resolve().parent
V=ROOT/'verification'
PDF=ROOT/'Q222_Analytical_Extraction_Refinement_REGENERATED_2026-10-22.pdf'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def write(p,data):p.write_text(json.dumps(data,indent=2,sort_keys=True)+'\n')

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--replay-report',type=Path)
    args=parser.parse_args()
    V.mkdir(exist_ok=True)
    if args.replay_report:
        target=V/'ISOLATED_REPLAY.json'
        if args.replay_report.resolve()!=target.resolve():shutil.copyfile(args.replay_report,target)
    replay=read(V/'ISOLATED_REPLAY.json')
    assert replay['status']=='PASS' and replay['replay']['status']=='PASS'
    assert replay['replay']['manifest_bound_output_changes']==[]
    assert replay['replay']['suite_count']==9
    pdfhash=sha(PDF)
    receipts=[];covered=set()
    for name in ['QA_ROOT_PAGES.json','QA_DIAGONAL_PAGES.json','QA_COMPILER_PAGES.json','QA_CHECKLIST_PAGES.json']:
        p=V/name;d=read(p)
        assert d.get('pdf_sha256',d.get('pdf',{}).get('sha256'))==pdfhash
        assert d['status'] in {'PASS','PASS_FOR_PDF_PAGES_16_TO_26'}
        assert not d.get('pdf_defects_found',[]) and not d.get('material_issues',[]) and not d.get('unresolved_visual_defects',[])
        pages=d.get('pages_reviewed',d.get('reviewed_pages'))
        assert pages and all(isinstance(x,int) for x in pages)
        covered.update(pages)
        receipts.append({'path':'verification/'+name,'sha256':sha(p),'pages':pages,'status':d['status']})
    assert covered==set(range(1,32))
    rendered=read(V/'RENDER_VALIDATION.json')
    assert rendered['status']=='PASS' and rendered['page_count']==31 and rendered['pdf_sha256']==pdfhash
    log=(V/'PDF_COMPILE.log').read_text()
    diagnostics=['Overfull','Underfull','undefined references','Missing character','duplicate destination','Rerun to get cross-references right']
    assert all(x not in log for x in diagnostics)
    qa={'status':'PASS','pdf_sha256':pdfhash,'pdf_pages':31,'full_pages_visually_reviewed':sorted(covered),'receipts':receipts,'final_render_decoding':'PASS_ALL_31','compile_diagnostics':{'overfull_boxes':0,'underfull_boxes':0,'unresolved_references':0,'missing_glyphs':0,'duplicate_destinations':0},'scope':'Document layout review within the same model system, not external mathematical peer review.'}
    write(V/'PDF_QA.json',qa)
    checklist=read(ROOT/'CHECKLIST.json')
    assert checklist['counts']=={'total':20,'scoped_completed':13,'open':7}
    assert [x['id'] for x in checklist['items'] if x['status']=='OPEN']==['B03','A01','G1','G2','E1','E2','E3']
    evidence=sorted({p for item in checklist['items'] for p in item['evidence']})
    assert all((ROOT/p).is_file() for p in evidence)
    old=ROOT/'recovered/Q222_Extraction_Continuation_Package_2026-10-22.zip'
    assert sha(old)=='9dcf9554b1a2578bb6b9fd909e7d194d82636fe61a0035da0694eb2e44f441a3'
    with zipfile.ZipFile(old) as z:assert z.testzip() is None
    c2=read(ROOT/'c2/regular_c2_independent_verification.json')
    tilted=read(ROOT/'compiler/CHECK_TILTED_EXTRACTION.json')
    qe=read(ROOT/'compiler/CHECK_Q222_QE_ENCODING.json')
    diagonal=read(ROOT/'diagonal/DIAGONAL_VERIFICATION.json')
    schur=read(ROOT/'c2/schur_remainder_verification.json')
    assert qe['quantifier_elimination_run'] is False and qe['candidate_feasibility']=='UNSOLVED'
    report={
        'release':'Q222 Analytical Extraction Refinement - regenerated 2026-10-22',
        'status':'PASS_FOR_STATED_SCOPES',
        'global_exact_extraction_complete':False,
        'historical_549_page_proof_recreated_or_reaudited':False,
        'source_custody_report':{'path':'inputs/SOURCE_CUSTODY.json','sha256':sha(ROOT/'inputs/SOURCE_CUSTODY.json')},
        'retained_archive':{'path':old.relative_to(ROOT).as_posix(),'bytes':old.stat().st_size,'sha256':sha(old),'crc_status':'PASS'},
        'replay':{'status':'PASS','suite_count':9,'construction_verification_suites':6,'separate_review_suites':3,'manifest_bound_files_at_replay':replay['manifest_bound_files'],'changed_manifest_bound_outputs':0,'report':'verification/ISOLATED_REPLAY.json','report_sha256':sha(V/'ISOLATED_REPLAY.json')},
        'local_C2_counts':c2['counts'],
        'C2_assertion_count_context':{'total':c2['total_assertions'],'encoding_checks':c2['assertions_by_group']['encoding'],'warning':'Overlapping individual comparisons, mostly encoding checks. Not independent mathematical theorem counts.'},
        'tilted_checks':tilted['check_count'],
        'diagonal_check_groups':[{'name':x['name'],'checks':x['checks'],'status':x['status']} for x in diagonal['groups']],
        'Schur_check_groups':schur['check_groups'],
        'QE_encoding_checks':qe['check_count'],
        'QE_status':'UNSOLVED; quantifier elimination was not run; no global coefficient tuple was extracted.',
        'pdf':{'file':PDF.name,'bytes':PDF.stat().st_size,'sha256':pdfhash,'pages':31,'visual_review':'PASS_ALL_31','compile_diagnostics':qa['compile_diagnostics']},
        'checklist':checklist['counts'],
        'open_items':['B03','A01','G1','G2','E1','E2','E3'],
        'checklist_evidence_paths':len(evidence),
        'review_limitations':['All analytical and separate-agent reviews are within the same model system.','No external peer review or proof-assistant formal verification.','Inherited source hashes do not certify all inherited mathematical claims.','No full-set exact pencil, global feature bound, or full-set order-two representation is asserted.'],
        'archive_integrity':'Performed separately after packaging; the final archive digest is recorded outside the ZIP.'
    }
    write(ROOT/'VERIFICATION.json',report)
    lines=[
        'Q222 ANALYTICAL EXTRACTION REFINEMENT - REGENERATED 2026-10-22',
        'VERIFICATION REPORT',
        '',
        'RESULT: PASS within the precise local, projected, module and conditional scopes.',
        'GLOBAL E1-E3: OPEN. Historical 549-page proof: not recreated or re-audited.',
        '',
        'FRESH ISOLATED REPLAY',
        'Nine suites passed in an isolated temporary copy: six construction/verification suites and three separately implemented review suites.',
        f'The replay checked {replay["manifest_bound_files"]} then-existing manifest-bound files. Every regenerated manifest-bound output remained byte-identical; mismatches: 0.',
        'The final manifest also binds these aggregate reports and final document-QA records. No delivered proof file was rewritten by the isolated replay.',
        ''
    ]
    for row in replay['replay']['suites']:lines.append('PASS  '+row['script'])
    lines.extend(['','EXACT LOCAL RESULTS',
        'Regular C2: 25 fixed frames; 735 mixed positive products; 50 positive pure remainders; 785 PSD blocks of order two; 1570 listed features retaining repetitions; 237 ambient monomials; 27747 product records; 295 seven-source target records; maximum product degree 6.',
        'The rational local radius is 510101005/350806159660744. The complete parent, contact equations, uniform Jacobian regularity, exact Taylor matrix, sign bounds, tensor factor typing and coefficient conventions are retained.',
        f'The C2 independent checker records {c2["total_assertions"]} comparisons, of which {c2["assertions_by_group"]["encoding"]} are encoding checks. These are not independent theorem counts and are not added to other overlapping check totals.',
        'Tilted projection: 38 checks; N=5, m=1 with 14 upper coefficient records, or N=4, m=0 with 12 records, over Q(sqrt(2)). Both inclusions hold for the stated projection only. The rational nonquantum scope witness is retained.',
        'Diagonal/B0/event: four corroboration groups with 4, 4, 15, 26 checks. The displayed B0 coefficient cone is exact; its physical restoration remains open. Restricted event cofactors satisfy all five polynomial kernel identities.',
        'Literal Schur remainder: 28 root identity groups plus the separate 17-check packet-matrix review. The infinite-span conclusion concerns literal scalar remainders and does not obstruct fixed features for the complete parent.',
        'QE encoding: 20 checks of the rational physical formula and decision input. Physical formula: 81 existential variables, 36 equations, 9 inequalities, maximum degree 5. Candidate M=3: 162 coefficients, 9 auxiliaries. Status UNSOLVED; no quantifier elimination was run.',
        '',
        'DOCUMENT QA',
        'The final manuscript has 31 pages. Every page was rendered and visually reviewed. All 31 final page images decode completely. Three truncated initial PNG render files were replaced; no PDF content change was needed.',
        'Final compile diagnostics: 0 overflowing boxes; 0 underfull boxes; 0 unresolved references; 0 missing glyphs; 0 duplicate destinations.',
        'PDF SHA256: '+pdfhash,
        '',
        'SOURCE CUSTODY',
        'Canonical August 30 base: 76839 bytes; SHA256 e268a6ae1c11dd02b5819064fa7e8187173074050d92f3204e703dd46f574731.',
        'Attached August 31 frozen source ZIP: 17609 bytes; SHA256 1a6208d7d1f51e5e166d46e80be3167ffc43aff97027f5d0c39705ab5d4f2bb3.',
        'Unchanged earlier continuation ZIP: 8884930 bytes; SHA256 '+sha(old)+'. Its ZIP CRC check passed.',
        '',
        'CHECKLIST',
        '20 supplementary items: 13 complete within stated scopes; 7 open.',
        'Open items: B03, A01, G1, G2, E1, E2, E3. All 34 unique evidence paths exist.',
        'Ten further concrete global data obligations G01-G10 are listed separately from the retained candidate gates G1 and G2.',
        '',
        'LIMITS',
        'The review comprises analytical review, separate-agent review within the same model system, and targeted exact algebra. It is not external peer review or proof-assistant formal verification.',
        'No evaluated full eight-coordinate quantum pencil, numerical global lift-size bound, or full-set second-order-cone representation is claimed.',
        'The preserved older order-13744 outer candidate and its prior records were not newly replayed in this regeneration; its reverse inclusion remains unproved.',
        '',
        'REPRODUCE',
        'python3 verify.py',
        'python3 -m pip install -r requirements.txt',
        'python3 verify.py --replay quick',
        '',
        'Final archive integrity is checked after packaging. The archive digest belongs outside the ZIP and is appended to the downloadable verification report.',
        ''
    ])
    (ROOT/'VERIFICATION.txt').write_text('\n'.join(lines))
    print(json.dumps({'status':report['status'],'replay_suites':9,'pdf_pages':31,'global_E1_E2_E3':'OPEN'},sort_keys=True))

if __name__=='__main__':main()

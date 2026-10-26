"""Bounded comparison of displayed equations in the scoped weak-C2 deltas."""
import collections,json,re
from pathlib import Path
HERE=Path(__file__).resolve().parent; ROOT=HERE.parent.parent
REF=ROOT/'recovery/refined_439/Q222_Refined_Proof_Source_2026-10-22'
manifest=json.loads((REF/'audit/weak_c2_scoping_manifest.json').read_text())
cross=json.loads((REF/'audit/component_crosswalk.json').read_text())
cross={x['component_id']:x for x in cross}
def equations(text):
    pat=r'\\\[.*?\\\]|\\begin\{(?:align\*?|equation\*?|gather\*?|multline\*?)\}.*?\\end\{(?:align\*?|equation\*?|gather\*?|multline\*?)\}'
    out=[]
    for match in re.finditer(pat,text,re.S):
        val=re.sub(r'\\(?:label|tag)\{[^}]*\}','',match.group())
        val=re.sub(r'\s+','',val)
        out.append(val)
    return collections.Counter(out)
report=[]
for row in manifest['fragments']:
    cid=row['component_id']; file=REF/cross[cid]['refined_file']
    new=file.read_text(encoding='utf-8')
    marker='\\label{component:'+cid+'}'
    start=new.index(marker)+len(marker)
    end=re.search(r'\\label\{component:p[0-9a-f]+\}',new[start:])
    body=new[start:start+end.start()] if end else new[start:]
    old=(ROOT/'recovery/full_r21/components/source'/f'{cid}.tex').read_text(encoding='utf-8')
    old_eq=equations(old); new_eq=equations(body)
    missing=list((old_eq-new_eq).elements())
    report.append({'component_id':cid,'refined_file':cross[cid]['refined_file'],'original_display_count':sum(old_eq.values()),'retained_original_display_count':sum(old_eq.values())-len(missing),'missing_displays':missing})
(HERE/'refined_display_comparison.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps({'compared_components':len(report),'missing_display_components':[{'component_id':x['component_id'],'missing':len(x['missing_displays']),'examples':[v[:130] for v in x['missing_displays'][:4]]} for x in report if x['missing_displays']]},indent=2))

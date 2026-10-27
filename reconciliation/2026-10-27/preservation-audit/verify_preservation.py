"""Audit Git-object preservation; this does not re-prove mathematical results."""
from pathlib import Path
import argparse
import json
import subprocess

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
DEFAULT_HEAD = 'ac8861aab22778ee06dbbf79c1784898e1e550f8'
BASES = ['1b293832ffeff46ebbee834b89554481d70f544e',
         '332483309586e0fec6e16fc1dd855915b648dca1']
PROTECTED = ('latex/', 'pdf/', 'extracted-lifts/',
             'complete-proof-materials/materials/', 'independent-review/checks/')


def git(*args):
    return subprocess.check_output(['git', *args], cwd=REPO)


def tree(ref):
    result = {}
    for row in git('ls-tree', '-r', '-z', ref).split(b'\0'):
        if row:
            meta, path = row.split(b'\t', 1)
            result[path.decode('utf-8')] = meta.decode().split()[2]
    return result


def document(ref, path):
    return json.loads(git('show', ref+':'+path))


parser = argparse.ArgumentParser()
parser.add_argument('--head', default=DEFAULT_HEAD)
parser.add_argument('--write-receipt', action='store_true')
args = parser.parse_args()
head = git('rev-parse', args.head).decode().strip()
current = tree(head)
current_checklist = document(head, 'complete-proof-materials/CHECKLIST.json')
reports = []
for base in BASES:
    subprocess.run(['git','merge-base','--is-ancestor',base,head],cwd=REPO,check=True)
    old = tree(base)
    changed = sorted(p for p in old.keys() & current.keys() if old[p] != current[p])
    removed = sorted(old.keys()-current.keys())
    old_list = document(base, 'complete-proof-materials/CHECKLIST.json')
    a = {r['id']:r for r in old_list['items']}
    b = {r['id']:r for r in current_checklist['items']}
    coverage = {key:old_list[key] == current_checklist[key]
                for key in ['components','c1_cells','c2_leaves']}
    protected = [p for p in old if p.startswith(PROTECTED)]
    report = {'base':base, 'base_tracked_files':len(old),
              'byte_identical_git_blobs':sum(old[p] == current.get(p) for p in old),
              'added_files':len(current.keys()-old.keys()), 'removed_files':removed,
              'changed_files':changed, 'protected_files':len(protected),
              'protected_changed_or_removed':[p for p in protected if old[p] != current.get(p)],
              'removed_checklist_ids':sorted(a.keys()-b.keys()),
              'changed_existing_checklist_items':[k for k in a if a[k] != b.get(k)],
              'coverage_records_identical':coverage}
    assert not removed and not report['protected_changed_or_removed']
    assert not report['removed_checklist_ids'] and all(coverage.values())
    if base == BASES[0]:
        assert not report['changed_existing_checklist_items']
    reports.append(report)

# The sole modified earlier mathematical candidate note has only added lines.
patch = git('diff','--no-ext-diff','--unified=0',BASES[1],head,'--',
            'independent-review/EXACT_LIFT_CANDIDATE.md').decode()
deleted_lines = [line for line in patch.splitlines() if line.startswith('-') and not line.startswith('---')]
assert not deleted_lines
old_m = document(BASES[0], 'complete-proof-materials/MANIFEST.json')
new_m = document(head, 'complete-proof-materials/MANIFEST.json')
a = {r['path']:r for r in old_m['files']}; b = {r['path']:r for r in new_m['files']}
assert not a.keys()-b.keys()
metadata_only = [p for p in a if a[p] != b[p] and
                 (a[p]['sha256'],a[p]['bytes']) == (b[p]['sha256'],b[p]['bytes'])]
result = {'scope':'Git-object and result-ledger preservation through the specified head; not a theorem audit',
          'head':head,'reports':reports,
          'earlier_candidate_note_has_no_deleted_lines':True,
          'latest_manifest_metadata_only_entries':metadata_only,
          'ok':True}
if args.write_receipt:
    (HERE/'result.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(result,indent=2))

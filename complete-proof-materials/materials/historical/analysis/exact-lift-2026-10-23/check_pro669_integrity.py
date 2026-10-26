from pathlib import Path
import hashlib, json, subprocess, zipfile

here = Path(__file__).resolve().parent
repo = here.parents[1] / 'github-q222-finite-sdp'
archive = Path('C:/Users/anubh/Downloads/Q222_Fixed_Family_Source_and_Matrices_2026-10-23.zip')
with zipfile.ZipFile(archive) as z:
    assert z.testzip() is None
    manifests = [n for n in z.namelist() if n.endswith('/MANIFEST.sha256')]
    assert len(manifests) == 1
    prefix = manifests[0][:-len('MANIFEST.sha256')]
    count = 0
    for line in z.read(manifests[0]).decode().splitlines():
        digest, name = line.split(maxsplit=1)
        assert hashlib.sha256(z.read(prefix + name)).hexdigest() == digest, name
        count += 1
    provenance = json.loads(z.read(prefix + 'root/INPUT_PROVENANCE.json'))
    checked = []
    for item in provenance['scientific_inputs']:
        data = subprocess.check_output(['git', 'show', provenance['commit'] + ':' + item['path']], cwd=repo)
        assert hashlib.sha256(data).hexdigest() == item['sha256'], item['path']
        checked.append(item['path'])
report = {'archive_sha256': hashlib.sha256(archive.read_bytes()).hexdigest(),
          'manifest_files_checked': count, 'input_commit': provenance['commit'],
          'git_input_hashes_checked': checked}
(here / 'pro669_integrity.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report, indent=2))

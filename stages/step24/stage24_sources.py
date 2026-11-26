from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parent

def stage(n):
    return ROOT.parent / f"step{n}"

def include(n):
    path = stage(n)
    if str(path) not in sys.path:
        sys.path.append(str(path))
    return path

def checked_files(n):
    from hashlib import sha256
    directory = stage(n)
    rows = (directory / 'files.sha256').read_text().splitlines()
    for row in rows:
        expected, name = row.split('  ', 1)
        if sha256((directory / name).read_bytes()).hexdigest() != expected:
            raise ValueError('Scientific source checksum mismatch: ' + name)
    return len(rows)

SHA = '94556ef2729ec2111154dcaf6a3e75891081cb6055d14d07bab0c09cfa3e8c2d'
def previous_modules():
    include(23)
    from stage23_sources import prepare
    old, model, data, info = prepare()
    return old, model, data, {'stage23_sha256': SHA, 'stage23_manifest_members': checked_files(23), 'inherited': info}

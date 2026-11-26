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

SHA = '18a9fcf5a7fd92555515a1e8eccf5d049cad5c242d7353f80d281558325f4a00'
def prepare():
    return include(16), checked_files(16)

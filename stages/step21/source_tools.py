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

SOURCE_SHA = 'b3006b1b524dc82ceb8ab2624cd638386ea7eac0cacb04aae7b56310f798df6b'
def prepare():
    checked_files(20)
    return stage(20)

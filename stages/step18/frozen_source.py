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

SHA17 = '8438f66de9c41d2845d5ddc712e2ef30e29f7f39de8bf4f3cdbaf8216be900aa'
def prepare():
    return include(17), include(16), checked_files(17), checked_files(16)

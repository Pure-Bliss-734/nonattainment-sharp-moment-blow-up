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

SOURCE_SHA = 'dab1de7256d730e4286dbf82cc38a0d7a136b53938e5a50918dbc2e9c10125e7'
PENCILS_SHA = '35f905894a2c61f6488bb0480c6c6551b229f63b5800a5bf92833dde277ba39e'
def prepare():
    return stage(21), stage(20), {f'stage{n}_members': checked_files(n) for n in (20,21)}

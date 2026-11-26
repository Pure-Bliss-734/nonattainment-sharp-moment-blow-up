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

SOURCE_SHA = 'f6415f4ddcdf874b3a738b9cff82e431c80fd7a83aa9e71ad45a0b83139f864d'
PENCILS_SHA = '35f905894a2c61f6488bb0480c6c6551b229f63b5800a5bf92833dde277ba39e'
def prepare():
    old=include(22)
    import bounded_model as bm
    return old, bm, bm.load_data(), {'stage22_sha256': SOURCE_SHA, 'stage22_manifest_members': checked_files(22), 'pencils_sha256': PENCILS_SHA}

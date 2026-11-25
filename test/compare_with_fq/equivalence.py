"""Runs fq lint (0.12.0) and dms_data_lint on 122 constructed FASTQ pairs and compares exit status and
first message (read errors by kind).
usage: python3 equivalence.py <fq binary> <dms_data_lint binary> <scratch dir>"""
import gzip, os, random, subprocess, sys, zlib

FQ, FL, D = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(D, exist_ok=True)
random.seed(1)

def rec(name, seq=None, qual=None, plus='+', n=20):
    seq = seq if seq is not None else ''.join(random.choice('ACGT') for _ in range(n))
    qual = qual if qual is not None else 'I' * len(seq)
    return f'{name}\n{seq}\n{plus}\n{qual}\n'

def pair(n=5, **kw):
    a = ''.join(rec(f'@r{i}/1') for i in range(n)); b = ''.join(rec(f'@r{i}/2') for i in range(n))
    return a, b

def gz(text, level=6):
    return gzip.compress(text.encode('latin-1') if isinstance(text, str) else text, compresslevel=level, mtime=0)

cases = []
def case(name, r1, r2, ext='.fq', raw=False):
    cases.append((name, r1, r2, ext, raw))

a, b = pair(50)
case('valid', a, b); case('valid gz', a, b, '.fq.gz')
case('valid gz level 1', gz(a, 1), gz(b, 1), '.fq.gz', True)
case('comments after space', ''.join(rec(f'@r{i} 1:N:0:ACGT') for i in range(5)), ''.join(rec(f'@r{i} 2:N:0:ACGT') for i in range(5)))
case('CRLF', a.replace('\n', '\r\n'), b.replace('\n', '\r\n'))
case('no final newline', a.rstrip('\n'), b.rstrip('\n'))
case('no final newline gz', a.rstrip('\n'), b.rstrip('\n'), '.fq.gz')
case('lowercase bases', a.replace('A', 'a'), b)
case('N bases', a.replace('G', 'N'), b)
case('all qualities 33-126', rec('@x/1', 'A' * 94, ''.join(chr(c) for c in range(33, 127))), rec('@x/2', 'A' * 94))
case('empty files', '', ''); case('empty gz members', '', '', '.fq.gz')
case('zero-byte gz files', b'', b'', '.fq.gz', True)
case('multi-member gz', gz(a[:400]) + gz(a[400:]), gz(b[:123]) + gz(b[123:]), '.fq.gz', True)
case('/1 vs /3 (equal after truncation)', rec('@r/1'), rec('@r/3'))
case('only one record', rec('@r/1'), rec('@r/2'))
big = 'ACGT' * (5 << 18)  # 5 MB sequence: longer than a chunk
case('record longer than a chunk', rec('@big/1', big) + a, rec('@big/2', big) + b, '.fq.gz')

def put(s, i, rep):  # replace record i of a pair text
    rs = [s[k:k + 0] for k in []]
    lines = s.split('\n'); lines[4 * i:4 * i + 4] = rep.rstrip('\n').split('\n'); return '\n'.join(lines)

case('S003 R1', put(a, 3, rec('r3/1')), b)
case('S003 R2', a, put(b, 3, rec('r3/2')))
case('S003 both, R1 reported', put(a, 3, rec('r3/1')), put(b, 3, rec('r3/2')))
case('S003 from truncation to nothing', put(a, 1, rec('/1')), put(b, 1, rec('/2')))
case('S004 empty sequence', put(a, 2, '@r2/1\n\n+\n\n'), b)
case('S004 empty plus line', put(a, 2, '@r2/1\nACGT\n\nIIII\n'), b)
case('S004 empty quality', a, put(b, 2, '@r2/2\nACGT\n+\n\n'))
case('S004 truncated final record', a + '@r50/1\nACGT\n', b + rec('@r50/2', 'ACGT'))
case('S004 name only at end', a + '@r50/1', b + '@r50/2')
case('trailing blank line', a + '\n', b)
case('trailing blank lines both', a + '\n\n', b + '\n\n')
case('S002 R in sequence', put(a, 4, rec('@r4/1', 'ACGRT')), b)
case('S002 dot', a, put(b, 4, rec('@r4/2', 'AC.GT')))
case('S002 non-ASCII byte', put(a, 4, rec('@r4/1', 'AC\xc3GT')), b)
case('S002 at position 40', put(a, 1, rec('@r1/1', 'A' * 39 + 'X' + 'A' * 20)), b)
case('S001 minus', put(a, 6, rec('@r6/1', plus='-')), b)
case('S001 plus with name', put(a, 6, rec('@r6/1', plus='+r6/1')), b)
case('S005 lengths', put(a, 7, rec('@r7/1', 'ACGT', 'III')), b)
case('S005 vs P001 in one record', put(a, 7, rec('@q7/1', 'ACGT', 'III')), b)
case('S006 space', put(a, 8, rec('@r8/1', 'ACGT', 'II I')), b)
case('S006 DEL', a, put(b, 8, rec('@r8/2', 'ACGT', 'II\x7fI')))
case('S006 tab', put(a, 8, rec('@r8/1', 'ACGT', 'I\tII')), b)
case('S002 R2 before S006 R1, same record', put(a, 8, rec('@r8/1', 'ACGT', 'II I')), put(b, 8, rec('@r8/2', 'ANGX')))
case('S006 R1 record 3 vs S002 R2 record 2', put(a, 3, rec('@r3/1', 'ACGT', 'II I')), put(b, 2, rec('@r2/2', 'ACXT')))
case('P001', put(a, 9, rec('@x9/1')), b)
case('P001 /1 x vs /2 x', rec('@r/1 x'), rec('@r/2 x'))
case('P001 differing comment', rec('@r 1:N'), rec('@r 1:Y'))
case('R1 shorter', a, b + rec('@r50/2'))
case('R2 shorter', a + rec('@r50/1'), b)
case('R2 shorter, gz', a + rec('@r50/1'), b, '.fq.gz')
dup = a + rec('@r7/1'); dupb = b + rec('@r7/2')
case('S007 duplicate', dup, dupb)
case('S007 duplicate gz', dup, dupb, '.fq.gz')
case('S007 duplicate after truncation', a + rec('@r7 extra'), b + rec('@r7 more'))
case('S007 two duplicates, first reported', a + rec('@r9/1') + rec('@r3/1'), b + rec('@r9/2') + rec('@r3/2'))
case('S007 vs S002 later', dup + rec('@z/1', 'AXGT'), dupb + rec('@z/2'))
# gzip damage
g = bytearray(gz(a * 1)); h = gz(b)
t = bytearray(g); t[len(t) // 2] ^= 0xFF
case('gz flipped byte', bytes(t), h, '.fq.gz', True)
t = bytearray(g); t[-8] ^= 0x01
case('gz bad CRC', bytes(t), h, '.fq.gz', True)
t = bytearray(g); t[-1] ^= 0x01
case('gz bad length', bytes(t), h, '.fq.gz', True)
case('gz truncated', bytes(g[:len(g) // 2]), h, '.fq.gz', True)
case('gz trailing garbage', bytes(g) + b'garbage', h, '.fq.gz', True)
case('gz trailing zeros', bytes(g) + b'\0' * 16, h, '.fq.gz', True)
case('plain text named .gz', a.encode(), h, '.fq.gz', True)
case('gzip named .fq (read as text)', bytes(g), b.encode(), '.fq', True)
case('R2 truncated gz', bytes(g), h[:len(h) - 30], '.fq.gz', True)
# a validation error before damage far later in the file
long_a = put(a, 1, rec('@r1/1', 'AXGT')) + ''.join(rec(f'@s{i}/1', n=100) for i in range(30000))
long_b = b + ''.join(rec(f'@s{i}/2', n=100) for i in range(30000))
t = bytearray(gz(long_a)); t[len(t) - 1000] ^= 0xFF
case('S002 at record 2, damage near the end', bytes(t), gz(long_b), '.fq.gz', True)
# random single-byte damage of plain text
for k in range(60):
    s = bytearray(a.encode()); pos = random.randrange(len(s)); s[pos] = random.choice(b'@+\n\r XN.\x00\x7fIAC')
    case(f'random edit {k} (byte {pos})', bytes(s), b.encode(), '.fq', True)

def write(path, data, ext, raw):
    with open(path, 'wb') as f:
        f.write(data if raw else (gz(data) if ext.endswith('.gz') else data.encode('latin-1')))

def run(cmd):
    p = subprocess.run(cmd, capture_output=True)
    err = [l for l in p.stderr.decode('utf-8', 'replace').splitlines() if ' INFO ' not in l and not l.startswith('dms_data_lint:')]
    return p.returncode, err

def key(rc, err):
    if rc == 0: return (0, '')
    first = err[0] if err else ''
    if first.startswith('Error: I/O error') or first.startswith('Error: could not'): return (rc, 'Error: I/O error')
    return (rc, first)

ok = 0; rows = []
for i, (name, r1, r2, ext, raw) in enumerate(cases):
    p1, p2 = f'{D}/c{i}_R1{ext}', f'{D}/c{i}_R2{ext}'
    write(p1, r1, ext, raw); write(p2, r2, ext, raw)
    f = key(*run([FQ, 'lint', p1, p2])); l = key(*run([FL, p1, p2]))
    same = f == l; ok += same
    rows.append((same, name, f, l))
    os.remove(p1); os.remove(p2)
for same, name, f, l in rows:
    if not same: print(f'DIFF  {name}\n   fq lint      : {f}\n   dms_data_lint: {l}')
print(f'{ok} of {len(cases)} cases identical (exit status and first error line; read errors compared by kind)')
for same, name, f, l in rows[:len(rows) - 60]:
    print(f"  {'=' if same else '!'} {name:45s} {f[1][:110] if f[0] else 'pass'}")

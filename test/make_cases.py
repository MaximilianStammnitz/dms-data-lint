"""Writes the small test cases in cases/ (each a pair <name>_1.fq / <name>_2.fq, or .fq.gz) from the
first 20 read pairs of the test data. The expected first message of each case, in cases/expected.tsv,
is the message of fq lint 0.12.0 on the same pair (read errors: only the kind, 'Error: I/O error').
usage: python3 make_cases.py   (from test/)"""
import gzip, os
os.makedirs('cases', exist_ok=True)
def recs(path, n=20):
    with gzip.open(path, 'rt') as f: lines = [next(f).rstrip('\n') for _ in range(4 * n)]
    return [lines[i:i + 4] for i in range(0, len(lines), 4)]
A, B = recs('GID1A_test_1.fastq.gz'), recs('GID1A_test_2.fastq.gz')
txt = lambda rs: ''.join('\n'.join(r) + '\n' for r in rs)
def edit(rs, i, k, new): rs = [list(r) for r in rs]; rs[i][k] = new; return rs
cases = {}
cases['valid'] = (txt(A), txt(B))
cases['crlf'] = (txt(A).replace('\n', '\r\n'), txt(B).replace('\n', '\r\n'))
cases['no_final_newline'] = (txt(A)[:-1], txt(B)[:-1])
cases['S003_missing_at'] = (txt(edit(A, 3, 0, A[3][0][1:])), txt(B))
cases['S004_empty_sequence'] = (txt(A), txt(edit(edit(B, 5, 1, ''), 5, 3, '')))
cases['S004_truncated_last_record'] = (txt(A) + A[0][0].replace(':', 'x:', 1) + '\nACGT\n', txt(B) + B[0][0].replace(':', 'x:', 1) + '\nACGT\n+\nIIII\n')
cases['S002_invalid_base'] = (txt(edit(A, 7, 1, A[7][1][:10] + 'R' + A[7][1][11:])), txt(B))
cases['S001_plus_line'] = (txt(edit(A, 8, 2, '-')), txt(B))
cases['S005_lengths'] = (txt(edit(A, 9, 3, A[9][3][:-1])), txt(B))
cases['S006_quality_char'] = (txt(A), txt(edit(B, 10, 3, B[10][3][:5] + ' ' + B[10][3][6:])))
cases['S002_R2_before_S006_R1'] = (txt(edit(A, 4, 3, ' ' + A[4][3][1:])), txt(edit(B, 4, 1, 'X' + B[4][1][1:])))
cases['P001_names'] = (txt(edit(A, 11, 0, '@another_read 1:N:0:1')), txt(B))
cases['r1_fewer_records'] = (txt(A[:19]), txt(B))
cases['r2_fewer_records'] = (txt(A), txt(B[:19]))
cases['S007_duplicate_name'] = (txt(A + [A[6]]), txt(B + [B[6]]))
gz = lambda s: gzip.compress(s.encode(), mtime=0)
binary = {}
binary['multi_member_gz'] = (gz(txt(A[:7])) + gz(txt(A[7:])), gz(txt(B)))
g = bytearray(gz(txt(A))); g[-8] ^= 1
binary['gz_bad_crc'] = (bytes(g), gz(txt(B)))
binary['gz_truncated'] = (gz(txt(A))[:-20], gz(txt(B)))
binary['gz_trailing_garbage'] = (gz(txt(A)) + b'junk', gz(txt(B)))
binary['not_gzip_named_gz'] = (txt(A).encode(), gz(txt(B)))
for name, (a, b) in cases.items():
    open(f'cases/{name}_1.fq', 'w').write(a); open(f'cases/{name}_2.fq', 'w').write(b)
for name, (a, b) in binary.items():
    open(f'cases/{name}_1.fq.gz', 'wb').write(a); open(f'cases/{name}_2.fq.gz', 'wb').write(b)

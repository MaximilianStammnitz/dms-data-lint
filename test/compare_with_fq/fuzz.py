"""Random damage: byte changes in gzip files and random edits of plain FASTQ text; fq lint vs
dms_data_lint. Prints, for plain and gzip files: [pairs, same pass/fail, same first message].
usage: python3 fuzz.py <fq binary> <dms_data_lint binary> <scratch dir> <number of pairs>"""
import gzip, os, random, subprocess, sys
FQ, FL, D, N = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]); os.makedirs(D, exist_ok=True); random.seed(7)
def rec(n, m): s=''.join(random.choice('ACGT') for _ in range(60)); return f'@r{n}:{random.randint(0,9)}/{m}\n{s}\n+\n{"I"*60}\n'
names=list(range(300)); A=''.join(rec(i,1) for i in names); B=''.join(rec(i,2) for i in names)
# make R2 names equal to R1 names
A_l=A.split('\n'); B_l=B.split('\n')
for i in range(0,len(A_l)-1,4): B_l[i]=A_l[i][:-1]+'2'
B='\n'.join(B_l)
def key(cmd):
    p=subprocess.run(cmd,capture_output=True); e=[l for l in p.stderr.decode('utf-8','replace').splitlines() if ' INFO ' not in l and not l.startswith('dms_data_lint:')]
    f=e[0] if e else ''
    return (p.returncode, 'Error: I/O error' if f.startswith('Error: I/O') else f) if p.returncode else (0,'')
diff=0; kinds={}; stat={'.fq':[0,0,0],'.fq.gz':[0,0,0]}
for k in range(N):
    gzmode = k % 2 == 0
    if gzmode:
        a=bytearray(gzip.compress(A.encode(),mtime=0)); b=gzip.compress(B.encode(),mtime=0)
        for _ in range(random.choice([1,1,2,3])): p=random.randrange(len(a)); a[p]=random.randrange(256)
        ext='.fq.gz'
    else:
        a=bytearray(A.encode()); b=B.encode()
        for _ in range(random.choice([1,2,4])):
            p=random.randrange(len(a)); op=random.random()
            if op<0.4: a[p]=random.choice(b'@+\n\r XN.\x00\x7f\tIAC/'); 
            elif op<0.7: del a[p]
            else: a[p:p]=bytes([random.choice(b'\n@+ /AI')])
        ext='.fq'
    p1,p2=f'{D}/f_R1{ext}',f'{D}/f_R2{ext}'; open(p1,'wb').write(a); open(p2,'wb').write(b)
    x=key([FQ,'lint',p1,p2]); y=key([FL,p1,p2])
    t=(x[1].split('] ')[0].split('[')[-1] if '[' in x[1] else x[1][:30]) or 'pass'; kinds[t]=kinds.get(t,0)+1
    st=stat[ext]; st[0]+=1; st[1]+=(x[0]==0)==(y[0]==0); st[2]+=x==y
    if x!=y:
        diff+=1
        if diff<=5 and ext=='.fq': print("DIFF", k, ext, "\n  fq      ", x, "\n  dms_data_lint", y)
print(stat); print(f'{N-diff} of {N} damaged pairs identical; outcomes:', dict(sorted(kinds.items(), key=lambda kv:-kv[1])))

#!/usr/bin/env python3
"""Black-box checks of opt-in grouping and unchanged exact-mode behavior."""
import csv
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

binary = Path(sys.argv[1]).resolve()
records = []

def read(name, exons, umi='AAAA', gene='G', cell='C', strand=0, source='S', contig='chr1', flag=0):
    cigar = ''.join((str(a-exons[i-1][1])+'N' if i else '') + str(b-a)+'M'
                    for i, (a,b) in enumerate(exons))
    length = sum(b-a for a,b in exons)
    return '\t'.join(map(str, [name, strand | flag, contig, exons[0][0]+1, 60, cigar,
                                '*', 0, 0, 'A'*length, 'I'*length,
                                'CB:Z:'+cell, 'UR:Z:'+umi, 'GX:Z:'+gene, 'SR:Z:'+source]))

def copies(prefix, count, exons, **kwargs):
    names = [prefix+str(i) for i in range(count)]
    records.extend(read(n, exons, **kwargs) for n in names)
    return set(names)

full = [(100,144),(250,300),(400,450)]
trunc = copies('full', 3, full, gene='trunc')
trunc |= copies('left',1,[(110,140)],gene='trunc')
trunc |= copies('right',1,[(420,440)],gene='trunc')
trunc |= copies('jitter',1,[(100,146),(248,302),(398,450)],gene='trunc')
long_chain = [(1000+200*i,1100+200*i) for i in range(14)]
long_group = copies('long',3,long_chain,gene='long')
long_group |= copies('prefix',1,long_chain[:3],gene='long')
long_group |= copies('suffix',1,long_chain[-3:],gene='long')
long_group |= copies('leftend',1,[(1010,1050)],gene='long')
long_group |= copies('rightend',1,[(3620,3650)],gene='long')
conflict_a = copies('a',3,full,gene='conflict')
conflict_b = copies('b',3,[full[0],full[-1]],gene='conflict')
bridge = copies('bridge',1,[(110,130)],gene='conflict')
retained = copies('retained',1,[(200,290)],gene='conflict')
# Counts from another structure must not irreversibly absorb a local UMI.
x = copies('x',100,full,gene='local',umi='AAAA')
y = copies('y_low',2,[full[0],full[-1]],gene='local',umi='AAAT')
y |= copies('y_high',3,[full[0],full[-1]],gene='local',umi='AATT')
# Keep an unsupported but genuinely conflicting isoform separate.
rare_a = copies('common',5,full,gene='rare')
rare_b = copies('rare',1,[full[0],full[-1]],gene='rare')
# A later conflicting anchor makes an earlier attachment ambiguous.
late_a = copies('late_a',3,[(500,600),(700,800),(900,1000)],gene='late')
late_c = copies('late_c',1,[(0,600),(700,750)],gene='late')
late_b = copies('late_b',3,[(500,600),(700,950)],gene='late')
# Observations outside the longest anchor can still give conflicting evidence.
outside_m = copies('outside_m',3,[(500,550),(650,700),(800,850)],gene='outside')
outside_a = copies('outside_a',3,[(100,150),(300,520)],gene='outside')
outside_b = copies('outside_b',3,[(100,520)],gene='outside')
# Deferred reads can contradict each other beyond both anchors' coverage.
central_p = [(500,550),(650,700),(800,850)]
central_q = [(500,550),(750,780),(800,850)]
deferred_p = copies('deferred_p',4,central_p,gene='deferred')
deferred_q = copies('deferred_q',4,central_q,gene='deferred')
deferred_c = copies('deferred_c',2,[(100,150),(300,520)],gene='deferred')
deferred_d = copies('deferred_d',2,[(100,520)],gene='deferred')
# Three mutually conflicting deferred geometries cannot fit into two anchors.
fallback_p = copies('fallback_p',4,central_p,gene='fallback')
fallback_q = copies('fallback_q',4,central_q,gene='fallback')
fallback_c = copies('fallback_c',2,[(100,150),(300,520)],gene='fallback')
fallback_e = copies('fallback_e',2,[(100,200),(350,520)],gene='fallback')
fallback_d = copies('fallback_d',2,[(100,520)],gene='fallback')
# A later witness can make two conflicting deferred nodes unique to one anchor.
unique_q = copies('unique_q',4,[(500,700),(800,850),(900,950)],gene='unique')
unique_p = copies('unique_p',4,central_p,gene='unique')
unique_c = copies('unique_c',2,[(100,150),(300,520)],gene='unique')
unique_d = copies('unique_d',2,[(100,200),(350,520)],gene='unique')
unique_e = copies('unique_e',2,[(100,600)],gene='unique')
# One junction cannot match two tiny introns merely because both are within tol.
short_a = copies('short_a',3,[(90,100),(102,200),(300,400),(500,520)],gene='short_n')
short_b = copies('short_b',3,[(90,95),(97,105),(107,115)],gene='short_n')
# Boundary dimensions must all remain effective.
scopes=[]
for label,kwargs in [('base',{}),('cell',{'cell':'C2'}),('strand',{'strand':16}),
                     ('source',{'source':'T'}),('contig',{'contig':'chr2'})]:
    scopes.append(copies('scope_'+label,2,full,gene='scope',**kwargs))
# Both derived tags and duplicate status are inherited by supplementary records.
records.append(read('full0',[(600,630)],gene='trunc',flag=2048))

header='@HD\tVN:1.6\tSO:unknown\n@SQ\tSN:chr1\tLN:100000\n@SQ\tSN:chr2\tLN:100000\n'

def partition(rows):
    groups=defaultdict(set)
    for row in rows: groups[row['molecule_id']].add(row['qname'])
    return set(map(frozenset,groups.values()))

def same_group(rows,names):
    ids={row['molecule_id'] for row in rows if row['qname'] in names}
    assert len(ids)==1,(names,ids)
    return next(iter(ids))

with tempfile.TemporaryDirectory(prefix='isoumi-compatible-') as td:
    root=Path(td)
    original=root/'input.sam'; original.write_text(header+'\n'.join(records)+'\n')
    reversed_input=root/'reverse.sam'; reversed_input.write_text(header+'\n'.join(reversed(records))+'\n')
    def run(label, inputs=(original,), flags=(), compatible=True, buckets=1, threads=1):
        prefix=root/label
        cmd=[str(binary),'--out',str(prefix),'--threads',str(threads),'--buckets',str(buckets),
             '--correction-method','directional','--source-tag','SR','--emit-tsv','--emit-explain','--mol-tag','MI']
        for inp in inputs: cmd += ['--bam',str(inp)]
        if compatible: cmd += ['--structure-mode','compatible']
        cmd += list(flags)
        subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=True)
        with prefix.with_suffix('.assignments.tsv').open() as f: rows=list(csv.DictReader(f,delimiter='\t'))
        with prefix.with_suffix('.corrections.tsv').open() as f: corrections=list(csv.DictReader(f,delimiter='\t'))
        return rows,corrections,prefix
    rows,corrections,prefix=run('compatible')
    assert len(rows)==len(records)-1
    same_group(rows,trunc); same_group(rows,long_group)
    ga=same_group(rows,conflict_a); gb=same_group(rows,conflict_b)
    assert ga != gb
    assert same_group(rows,bridge) in {ga,gb}
    assert next(r for r in rows if r['qname']=='bridge0')['structure_status']=='ambiguous'
    assert same_group(rows,retained) not in {ga,gb}
    assert next(r for r in rows if r['qname']=='retained0')['structure_status']=='unsupported'
    assert same_group(rows,x) != same_group(rows,y)
    ratio_rows,_,_=run('ratio',flags=['--correction-method','ratio'])
    assert len({r['molecule_id'] for r in ratio_rows if r['qname'] in y})==2
    assert same_group(ratio_rows,x) not in {r['molecule_id'] for r in ratio_rows if r['qname'] in y}
    assert any(r['raw_umi']=='AAAT' and r['corr_umi']=='AATT' and r['raw_count']=='2'
               and r['seed_count']=='3' for r in corrections)
    assert same_group(rows,rare_a) != same_group(rows,rare_b)
    assert same_group(rows,late_a) != same_group(rows,late_b)
    assert next(r for r in rows if r['qname']=='late_c0')['structure_status']=='ambiguous'
    assert same_group(rows,outside_a) != same_group(rows,outside_b)
    assert same_group(rows,deferred_p) != same_group(rows,deferred_q)
    assert same_group(rows,deferred_c) != same_group(rows,deferred_d)
    assert all(r['structure_status']=='ambiguous' for r in rows
               if r['qname'] in deferred_c | deferred_d)
    assert len({same_group(rows,names) for names in [fallback_c,fallback_e,fallback_d]})==3
    assert same_group(rows,fallback_d) not in {same_group(rows,fallback_p),same_group(rows,fallback_q)}
    assert all(r['structure_status']=='unsupported' for r in rows if r['qname'] in fallback_d)
    assert same_group(rows,unique_p) != same_group(rows,unique_q)
    assert same_group(rows,unique_c) != same_group(rows,unique_d)
    assert same_group(rows,unique_e) == same_group(rows,unique_q)
    assert all(r['structure_status']=='unsupported' for r in rows if r['qname'] in unique_d)
    assert same_group(rows,short_a) != same_group(rows,short_b)
    assert all(r['structure_status']=='compatible' for r in rows if r['qname'] in short_a | short_b)
    assert len({same_group(rows,names) for names in scopes})==len(scopes)
    reversed_rows,_,_=run('reversed',(reversed_input,),buckets=7,threads=2)
    assert partition(reversed_rows)==partition(rows)
    exact,_,_=run('exact',compatible=False)
    explicit,_,_=run('explicit',compatible=False,flags=['--structure-mode','exact'])
    assert partition(exact)==partition(explicit)
    assert 'structure_status' not in exact[0]
    assert len({r['molecule_id'] for r in exact if r['qname'] in trunc})>1
    # Input isolation remains independent of identical CB/GX/UMI/geometry.
    left=root/'one.sam'; right=root/'two.sam'
    left.write_text(header+read('one',full)+'\n'); right.write_text(header+read('two',full)+'\n')
    shared,_,_=run('shared',(left,right)); isolated,_,_=run('isolated',(left,right),flags=['--isolate-inputs'])
    assert len(partition(shared))==1 and len(partition(isolated))==2
    for flags in [('--no-structure',),('--end-bin','100')]:
        cmd=[str(binary),'--bam',str(original),'--out',str(root/'invalid'),
             '--structure-mode','compatible',*flags]
        assert subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE).returncode != 0
    samtools=shutil.which('samtools')
    if samtools:
        raw=subprocess.check_output([samtools,'view',str(prefix)+'.dedup.bam'],text=True)
        alignments=[]
        for line in raw.splitlines():
            fields=line.split('\t'); tags={x[:2]:x[5:] for x in fields[11:]}
            alignments.append((fields[0],int(fields[1]),tags))
        assert len(alignments)==len(records)
        primary={n:(flag,tags) for n,flag,tags in alignments if not (flag & 2304)}
        supplements=[(n,flag,tags) for n,flag,tags in alignments if flag & 2048]
        assert len(supplements)==1
        n,flag,tags=supplements[0]
        for tag in ['UB','MI','DA']: assert tags[tag]==primary[n][1][tag]
        assert bool(flag & 1024)==bool(primary[n][0] & 1024)
        for n in y: assert primary[n][1]['UB']=='AATT'
        # Every final corrected UMI has a raw-UMI observation in that final group.
        by_mi=defaultdict(list)
        for n,(flag,tags) in primary.items(): by_mi[tags['MI']].append((n,flag,tags))
        by_q={r['qname']:r for r in rows}
        for group in by_mi.values():
            assert len({tags['UB'] for n,flag,tags in group})==1
            assert group[0][2]['UB'] in {tags['UR'] for n,flag,tags in group}
            assert sum(tags['DA']=='0' for n,flag,tags in group)==1
            for n,flag,tags in group:
                assert tags['MI']==by_q[n]['molecule_id']
                assert tags['DA']==by_q[n]['dup']==str(int(bool(flag & 1024)))
    else:
        print('SKIP: samtools unavailable; compatible BAM-tag inheritance subchecks not run',file=sys.stderr)
print('PASS: compatible structure grouping, local UMI counts, scopes, ambiguity, and exact-mode parity')

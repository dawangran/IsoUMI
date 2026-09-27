#!/usr/bin/env python3
"""Mixed alignment directions must not split a gene-assigned RNA molecule."""
import csv
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path


binary = Path(sys.argv[1]).resolve()
samtools = shutil.which('samtools')
header = ('@HD\tVN:1.6\tSO:unknown\n'
          '@SQ\tSN:chr1\tLN:100000\n@SQ\tSN:chr2\tLN:100000\n')
full = [(100, 150), (250, 300), (400, 450)]
structures = {'exact': (), 'compatible': ('--structure-mode', 'compatible'),
              'no_structure': ('--no-structure',)}


def read(name, flag=0, gene='G', gene_type='Z', cell='C', umi='AAAA',
         exons=full, contig='chr1', extra=()):
    length = sum(end-start for start, end in exons)
    cigar = ''.join((str(start-exons[i-1][1])+'N' if i else '') +
                    str(end-start)+'M' for i, (start, end) in enumerate(exons))
    seq = ('ACGTTGCAAA' * ((length+9)//10))[:length]
    tags = ['CB:Z:'+cell, 'UR:Z:'+umi]
    if gene is not None:
        tags.append('GX:'+gene_type+':'+gene)
    if flag & 4:
        contig, pos, cigar, mapq = '*', 0, '*', 0
    else:
        pos, mapq = exons[0][0]+1, 60
    return '\t'.join(map(str, [name, flag, contig, pos, mapq, cigar, '*', 0, 0,
                                seq, 'I'*length, *tags, *extra]))


def pair(prefix, **kwargs):
    return [read(prefix+'_f', 0, **kwargs), read(prefix+'_r', 16, **kwargs)]


def ids(rows, names):
    return {row['molecule_id'] for row in rows if row['qname'] in names}


def names(records):
    return {record.split('\t', 1)[0] for record in records}


def table(path):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


with tempfile.TemporaryDirectory(prefix='isoumi-strand-') as td:
    root = Path(td)

    def run(label, records, flags=(), buckets=1, threads=1):
        source = root/(label+'.sam')
        source.write_text(header+'\n'.join(records)+'\n')
        prefix = root/label
        command = [str(binary), '--bam', str(source), '--out', str(prefix),
                   '--threads', str(threads), '--buckets', str(buckets),
                   '--emit-tsv', '--mol-tag', 'MI', *flags]
        result = subprocess.run(command, capture_output=True, text=True)
        assert result.returncode == 0, (label, result.stderr)
        rows = table(Path(str(prefix)+'.assignments.tsv'))
        molecules = table(Path(str(prefix)+'.molecules.tsv'))
        alignments = []
        if samtools:
            raw = subprocess.check_output([samtools, 'view', str(prefix)+'.dedup.bam'], text=True)
            originals = {(fields[0], int(fields[1])): fields
                         for fields in (record.split('\t') for record in records)}
            for line in raw.splitlines():
                fields = line.split('\t')
                flag = int(fields[1])
                original = originals[(fields[0], flag & ~1024)]
                # No orientation, coordinates, CIGAR, sequence, or quality rewrite.
                assert fields[2:11] == original[2:11], (label, fields[0])
                tags = {tag[:2]: tag[5:] for tag in fields[11:]}
                original_tags = {tag[:2]: tag[5:] for tag in original[11:]}
                assert tags['UR'] == original_tags['UR']
                alignments.append((fields[0], flag, tags))
            assert len(alignments) == len(records)
            by_qname = {row['qname']: row for row in rows}
            for qname, flag, tags in alignments:
                if flag & (4 | 256 | 2048):
                    continue
                assert tags['MI'] == by_qname[qname]['molecule_id']
                assert tags['DA'] == by_qname[qname]['dup'] == str(int(bool(flag & 1024)))
        return rows, molecules, alignments

    # Four forward, four reverse, four mixed molecules; UMIs are >=4 edits apart.
    truth = []
    truth_names = []
    for molecule in range(12):
        umi = ''.join('ACGT'[(molecule//(4**digit)) % 4]*4 for digit in range(3))
        count = 20 if molecule < 8 else 25
        records = []
        for index in range(count):
            flag = 0 if molecule < 4 else 16 if molecule < 8 else (index % 2)*16
            records.append(read('m%02d_r%02d' % (molecule, index), flag, umi=umi))
        truth.extend(records)
        truth_names.append(names(records))
    assert len(truth) == 260

    for structure, structure_flags in structures.items():
        for method in ('ratio', 'directional'):
            baseline = None
            for strand in ('default', 'auto', 'ignore', 'alignment'):
                flags = (*structure_flags, '--correction-method', method)
                if strand != 'default':
                    flags += ('--strand-mode', strand)
                label = '_'.join((structure, method, strand))
                rows, molecules, alignments = run(label, truth, flags)
                expected = 16 if strand == 'alignment' else 12
                assert len(rows) == 260 and len(molecules) == expected, label
                assert sum(int(m['count']) for m in molecules) == 260
                all_ids = set()
                for index, group_names in enumerate(truth_names):
                    group_ids = ids(rows, group_names)
                    assert len(group_ids) == (2 if strand == 'alignment' and index >= 8 else 1), label
                    assert not all_ids & group_ids
                    all_ids |= group_ids
                    counts = Counter(row['molecule_id'] for row in rows if row['qname'] in group_names)
                    expected_counts = [12, 13] if len(group_ids) == 2 else [len(group_names)]
                    assert sorted(counts.values()) == expected_counts
                    for molecule_id in group_ids:
                        assert sum(row['dup'] == '0' for row in rows
                                   if row['molecule_id'] == molecule_id) == 1
                assert sorted(int(m['count']) for m in molecules) == sorted(Counter(
                    row['molecule_id'] for row in rows).values())
                if strand != 'alignment':
                    assert all('|STR=.|' in m['key'] for m in molecules)
                normalized = {r['qname']: (r['molecule_id'], r['dup']) for r in rows}
                if strand == 'default':
                    baseline = normalized
                elif strand in ('auto', 'ignore'):
                    assert normalized == baseline
            reversed_rows, _, _ = run(structure+'_'+method+'_reversed', list(reversed(truth)),
                                      (*structure_flags, '--correction-method', method),
                                      buckets=7, threads=2)
            # Equal-quality representatives use input order; molecule IDs do not.
            assert {r['qname']: r['molecule_id'] for r in reversed_rows} == {
                qname: value[0] for qname, value in baseline.items()}
            parallel_rows, _, _ = run(structure+'_'+method+'_parallel', truth,
                                      (*structure_flags, '--correction-method', method),
                                      buckets=7, threads=2)
            assert {r['qname']: (r['molecule_id'], r['dup']) for r in parallel_rows} == baseline

    # Correct a low-frequency one-edit UMI after pooling both orientations.
    erroneous = [read('seed_%02d' % i, (i % 2)*16, umi='AAAA') for i in range(20)]
    erroneous += [read('error_%02d' % i, i*16, umi='AAAT') for i in range(2)]
    for structure, flags in structures.items():
        for method in ('ratio', 'directional'):
            rows, molecules, alignments = run('correction_'+structure+'_'+method, erroneous,
                                             (*flags, '--correction-method', method))
            assert len(rows) == 22 and len(molecules) == 1
            assert molecules[0]['umi_corr'] == 'AAAA' and molecules[0]['count'] == '22'
            assert len(ids(rows, names(erroneous))) == 1
            assert sum(row['dup'] == '0' for row in rows) == 1
            assert all(tags['UB'] == 'AAAA' for _, _, tags in alignments)

    # Automatic fallback must inspect the selected tag's actual string value.
    fallback = []
    fallback_names = {}
    for label, kwargs in [('present', {}), ('missing', {'gene': None}),
                          ('empty', {'gene': ''}), ('numeric', {'gene': '7', 'gene_type': 'i'}),
                          ('hex', {'gene': 'ABCD', 'gene_type': 'H'}),
                          ('literal_na', {'gene': 'NA'})]:
        records = pair(label, cell=label, **kwargs)
        fallback.extend(records)
        fallback_names[label] = names(records)
    for structure, flags in structures.items():
        for strand in ('auto', 'alignment', 'ignore'):
            rows, _, _ = run('fallback_'+structure+'_'+strand, fallback,
                             (*flags, '--strand-mode', strand))
            for label, group_names in fallback_names.items():
                expected = 1 if strand == 'ignore' or (strand == 'auto' and label in ('present', 'literal_na')) else 2
                assert len(ids(rows, group_names)) == expected, (structure, strand, label)
        for strand, expected in [('auto', 2), ('ignore', 1), ('alignment', 2)]:
            rows, _, _ = run('no_gene_'+structure+'_'+strand, pair('no_gene'),
                             (*flags, '--no-gene', '--strand-mode', strand))
            assert len(ids(rows, names(pair('no_gene')))) == expected

        # GX differs but GN agrees: --gene-tag must control auto mode and scope.
        custom = [read('gn_f', 0, gene='X', extra=['GN:Z:G']),
                  read('gn_r', 16, gene='Y', extra=['GN:Z:G']),
                  *pair('gx_only', cell='GX_only')]
        rows, _, _ = run('custom_gene_'+structure, custom, (*flags, '--gene-tag', 'GN'))
        assert len(ids(rows, {'gn_f', 'gn_r'})) == 1
        assert len(ids(rows, {'gx_only_f', 'gx_only_r'})) == 2

        # Ignoring orientation must preserve independent biological boundaries.
        scopes = {}
        for label, kwargs in [('base', {}), ('cell', {'cell': 'C2'}), ('gene', {'gene': 'G2'}),
                              ('umi', {'umi': 'TTTT'}), ('contig', {'contig': 'chr2'}),
                              ('structure', {'exons': [full[0], full[-1]]})]:
            scopes[label] = pair('scope_'+label, **kwargs)
        rows, _, _ = run('boundaries_'+structure, sum(scopes.values(), []), flags)
        group_ids = {label: ids(rows, names(records)) for label, records in scopes.items()}
        assert all(len(value) == 1 for value in group_ids.values())
        expected = 5 if structure == 'no_structure' else 6
        assert len(set.union(*group_ids.values())) == expected
        assert (group_ids['base'] == group_ids['structure']) == (structure == 'no_structure')

    # End bins use genomic left/right coordinates when orientation is ignored.
    endpoints = [*pair('ends'), *pair('left_shift', exons=[(40, 150), *full[1:]]),
                 *pair('right_shift', exons=[*full[:-1], (400, 510)])]
    for strand in ('auto', 'ignore', 'alignment'):
        rows, molecules, _ = run('ends_'+strand, endpoints,
                                 ('--strand-mode', strand, '--end-bin', '50'))
        assert len(molecules) == (6 if strand == 'alignment' else 3)
        if strand == 'alignment':
            assert all('|E5=' in row['key'] and '|E3=' in row['key'] for row in molecules)
        else:
            assert len(ids(rows, {'ends_f', 'ends_r'})) == 1
            assert any('|EL=100|ER=450' in row['key'] for row in molecules)
            assert any('|EL=0|ER=450' in row['key'] for row in molecules)
            assert any('|EL=100|ER=500' in row['key'] for row in molecules)
            assert all('|E5=' not in row['key'] for row in molecules)

    # Non-primary records inherit their primary read's final cross-strand status.
    inherited = [*pair('inherit'), read('inherit_f', 2048, exons=[(700, 720)]),
                 read('inherit_r', 256, exons=[(900, 920)]),
                 read('inherit_r', 2048 | 16, exons=[(1100, 1120)]),
                 read('unmapped', 4 | 16, extra=['MI:Z:stale', 'DA:i:1'])]
    for structure, flags in structures.items():
        rows, molecules, alignments = run('inherit_'+structure, inherited, flags, buckets=3)
        assert len(rows) == 2 and len(molecules) == 1 and molecules[0]['count'] == '2'
        if samtools:
            primary = {qname: (flag, tags) for qname, flag, tags in alignments
                       if not flag & (4 | 256 | 2048)}
            for qname, flag, tags in alignments:
                if flag & 4:
                    assert 'MI' not in tags and tags['DA'] == '0' and tags['UB'] == tags['UR']
                    assert not flag & 1024
                elif flag & (256 | 2048):
                    pflag, ptags = primary[qname]
                    assert all(tags[tag] == ptags[tag] for tag in ('MI', 'UB', 'DA'))
                    assert bool(flag & 1024) == bool(pflag & 1024)

    invalid = subprocess.run([str(binary), '--strand-mode', 'invalid'], capture_output=True, text=True)
    assert invalid.returncode != 0 and '--strand-mode' in invalid.stderr

if not samtools:
    print('SKIP: samtools unavailable; strand BAM preservation and inheritance subchecks not run', file=sys.stderr)
print('PASS: mixed-strand counting, strand modes, gene fallback, boundaries, end bins, and determinism')

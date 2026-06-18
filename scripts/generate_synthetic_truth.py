#!/usr/bin/env python3
"""Generate a deterministic synthetic truth dataset for IsoUMI benchmarks.

The dataset is intentionally small but covers the core Application Note claims:

- low-support UMI errors should collapse within the same transcript context;
- identical UMIs in different splice-junction chains should not collapse;
- identical UMIs in different cells or chromosomes should remain separate;
- non-spliced reads use locus-bin grouping;
- reads missing cell barcodes pass through without correction reports.
"""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path


SAM_HEADER = (
    "@HD\tVN:1.6\tSO:coordinate\n"
    "@SQ\tSN:chr1\tLN:100000\n"
    "@SQ\tSN:chr2\tLN:100000\n"
)


@dataclass(frozen=True)
class Template:
    scenario: str
    true_molecule: str
    cell: str | None
    gene: str
    chrom: str
    pos: int
    strand: str
    cigar: str
    true_umi: str
    raw_umi: str
    umi_qual: str
    expected_corr_umi: str
    group_class: str
    naive_group: str
    should_report: bool = True


def query_length(cigar: str) -> int:
    length = 0
    for n_str, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        n = int(n_str)
        if op in {"M", "I", "S", "=", "X"}:
            length += n
    if length <= 0:
        raise ValueError(f"CIGAR does not consume query bases: {cigar}")
    return length


def repeated_sequence(length: int) -> str:
    motif = "ACGT"
    return (motif * ((length // len(motif)) + 1))[:length]


def sam_record(qname: str, tmpl: Template) -> str:
    qlen = query_length(tmpl.cigar)
    seq = repeated_sequence(qlen)
    qual = "F" * qlen
    flag = 16 if tmpl.strand == "-" else 0
    tags = []
    if tmpl.cell is not None:
        tags.append(f"CB:Z:{tmpl.cell}")
    tags.extend(
        [
            f"UR:Z:{tmpl.raw_umi}",
            f"UY:Z:{tmpl.umi_qual}",
            f"GX:Z:{tmpl.gene}",
            f"ZT:Z:{tmpl.true_molecule}",
            f"ZU:Z:{tmpl.true_umi}",
            f"ZS:Z:{tmpl.scenario}",
        ]
    )
    fields = [
        qname,
        str(flag),
        tmpl.chrom,
        str(tmpl.pos),
        "60",
        tmpl.cigar,
        "*",
        "0",
        "0",
        seq,
        qual,
        *tags,
    ]
    return "\t".join(fields)


def truth_row(qname: str, tmpl: Template) -> str:
    cell = tmpl.cell if tmpl.cell is not None else "NA"
    fields = [
        qname,
        tmpl.scenario,
        tmpl.true_molecule,
        cell,
        tmpl.gene,
        tmpl.chrom,
        str(tmpl.pos),
        tmpl.strand,
        tmpl.cigar,
        tmpl.true_umi,
        tmpl.raw_umi,
        tmpl.umi_qual,
        tmpl.expected_corr_umi,
        tmpl.group_class,
        tmpl.naive_group,
        "1" if tmpl.should_report else "0",
    ]
    return "\t".join(fields)


def add_replicates(records: list[Template], tmpl: Template, n: int) -> None:
    records.extend([tmpl] * n)


def build_templates() -> list[Template]:
    records: list[Template] = []

    # Same cell/gene/chromosome/junction: one low-quality 1-mismatch UMI should
    # collapse into the high-support AAAA molecule.
    add_replicates(
        records,
        Template(
            scenario="umi_error_same_junction",
            true_molecule="MOL_JUNC_A",
            cell="CELL_A",
            gene="GENE1",
            chrom="chr1",
            pos=101,
            strand="+",
            cigar="50M100N50M",
            true_umi="AAAA",
            raw_umi="AAAA",
            umi_qual="IIII",
            expected_corr_umi="AAAA",
            group_class="chr1_plus_junction_A",
            naive_group="CELL_A|GENE1|AAAA",
        ),
        5,
    )
    add_replicates(
        records,
        Template(
            scenario="umi_error_same_junction",
            true_molecule="MOL_JUNC_A",
            cell="CELL_A",
            gene="GENE1",
            chrom="chr1",
            pos=101,
            strand="+",
            cigar="50M100N50M",
            true_umi="AAAA",
            raw_umi="AAAT",
            umi_qual="!!!!",
            expected_corr_umi="AAAA",
            group_class="chr1_plus_junction_A",
            naive_group="CELL_A|GENE1|AAAT",
        ),
        1,
    )

    # Same cell/gene/raw UMI but different junction chains: a naive
    # cell+gene+UMI baseline over-collapses these, while IsoUMI should not.
    add_replicates(
        records,
        Template(
            scenario="same_umi_different_junction",
            true_molecule="MOL_JUNC_B",
            cell="CELL_A",
            gene="GENE1",
            chrom="chr1",
            pos=101,
            strand="+",
            cigar="50M200N50M",
            true_umi="CCCC",
            raw_umi="CCCC",
            umi_qual="IIII",
            expected_corr_umi="CCCC",
            group_class="chr1_plus_junction_B",
            naive_group="CELL_A|GENE1|CCCC",
        ),
        3,
    )
    add_replicates(
        records,
        Template(
            scenario="same_umi_different_junction",
            true_molecule="MOL_JUNC_C",
            cell="CELL_A",
            gene="GENE1",
            chrom="chr1",
            pos=101,
            strand="+",
            cigar="50M300N50M",
            true_umi="CCCC",
            raw_umi="CCCC",
            umi_qual="IIII",
            expected_corr_umi="CCCC",
            group_class="chr1_plus_junction_C",
            naive_group="CELL_A|GENE1|CCCC",
        ),
        3,
    )

    # Same UMI in a different cell and chromosome should remain separate.
    add_replicates(
        records,
        Template(
            scenario="same_umi_different_cell",
            true_molecule="MOL_CELL_B",
            cell="CELL_B",
            gene="GENE1",
            chrom="chr1",
            pos=101,
            strand="+",
            cigar="50M100N50M",
            true_umi="AAAA",
            raw_umi="AAAA",
            umi_qual="IIII",
            expected_corr_umi="AAAA",
            group_class="cell_b_chr1_plus_junction_A",
            naive_group="CELL_B|GENE1|AAAA",
        ),
        2,
    )
    add_replicates(
        records,
        Template(
            scenario="same_umi_different_chromosome",
            true_molecule="MOL_CHR2",
            cell="CELL_A",
            gene="GENE1",
            chrom="chr2",
            pos=101,
            strand="+",
            cigar="50M100N50M",
            true_umi="AAAA",
            raw_umi="AAAA",
            umi_qual="IIII",
            expected_corr_umi="AAAA",
            group_class="chr2_plus_junction_A",
            naive_group="CELL_A|GENE1|AAAA",
        ),
        2,
    )

    # Non-spliced reads exercise the locus-bin path.
    add_replicates(
        records,
        Template(
            scenario="umi_error_locus_bin",
            true_molecule="MOL_LOCUS_A",
            cell="CELL_A",
            gene="GENE2",
            chrom="chr1",
            pos=1001,
            strand="+",
            cigar="100M",
            true_umi="TTTT",
            raw_umi="TTTT",
            umi_qual="IIII",
            expected_corr_umi="TTTT",
            group_class="chr1_plus_locus_bin_A",
            naive_group="CELL_A|GENE2|TTTT",
        ),
        4,
    )
    add_replicates(
        records,
        Template(
            scenario="umi_error_locus_bin",
            true_molecule="MOL_LOCUS_A",
            cell="CELL_A",
            gene="GENE2",
            chrom="chr1",
            pos=1001,
            strand="+",
            cigar="100M",
            true_umi="TTTT",
            raw_umi="TTTA",
            umi_qual="!!!!",
            expected_corr_umi="TTTT",
            group_class="chr1_plus_locus_bin_A",
            naive_group="CELL_A|GENE2|TTTA",
        ),
        1,
    )

    # Missing cell barcodes should pass through and not appear in correction or
    # assignment reports.
    add_replicates(
        records,
        Template(
            scenario="missing_cell_passthrough",
            true_molecule="MOL_NO_CELL",
            cell=None,
            gene="GENE1",
            chrom="chr1",
            pos=101,
            strand="+",
            cigar="50M100N50M",
            true_umi="GGGG",
            raw_umi="GGGA",
            umi_qual="!!!!",
            expected_corr_umi="GGGA",
            group_class="missing_cell",
            naive_group="NA|GENE1|GGGA",
            should_report=False,
        ),
        1,
    )

    return records


def write_dataset(out_dir: Path, prefix: str) -> tuple[Path, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    sam_path = out_dir / f"{prefix}.sam"
    truth_path = out_dir / f"{prefix}.truth.tsv"
    templates = build_templates()

    with sam_path.open("w", encoding="utf-8") as sam:
        sam.write(SAM_HEADER)
        for idx, tmpl in enumerate(templates, start=1):
            qname = f"sim{idx:04d}_{tmpl.scenario}"
            sam.write(sam_record(qname, tmpl))
            sam.write("\n")

    with truth_path.open("w", encoding="utf-8") as truth:
        truth.write(
            "qname\tscenario\ttrue_molecule\tcell\tgene\tchrom\tpos\tstrand\t"
            "cigar\ttrue_umi\traw_umi\tumi_qual\texpected_corr_umi\t"
            "group_class\tnaive_group\tshould_report\n"
        )
        for idx, tmpl in enumerate(templates, start=1):
            qname = f"sim{idx:04d}_{tmpl.scenario}"
            truth.write(truth_row(qname, tmpl))
            truth.write("\n")

    return sam_path, truth_path


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Generate a deterministic synthetic truth SAM for IsoUMI."
    )
    parser.add_argument(
        "--out-dir",
        default="benchmarks/output/synthetic_truth",
        help="Output directory (default: benchmarks/output/synthetic_truth).",
    )
    parser.add_argument(
        "--prefix",
        default="synthetic_truth",
        help="Output filename prefix (default: synthetic_truth).",
    )
    args = parser.parse_args()

    sam_path, truth_path = write_dataset(Path(args.out_dir), args.prefix)
    print(f"Wrote SAM: {sam_path}")
    print(f"Wrote truth: {truth_path}")
    print("Suggested IsoUMI command:")
    print(
        f"  src/isoumi --bam {sam_path} --out {Path(args.out_dir) / 'isoumi'} "
        "--threads 1 --buckets 8 --ham 1 --ratio 0.60 --emit-tsv --emit-explain"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

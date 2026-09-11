#!/usr/bin/env python3
"""Collate kir-mapper's nested donor output into stable, flat tables."""

from __future__ import annotations

import argparse
import csv
import gzip
import tarfile
from pathlib import Path


CALL_COLUMNS = ["Sample", "Copy_number", "Calls", "Ratio", "Missings"]
REPORT_COLUMNS = [
    "Sample",
    "Copy_number",
    "Chr",
    "Allele_A",
    "Allele_B",
    "Allele_C",
    "Allele_D",
    "Tested_genotypes",
    "Valid_genotypes",
    "Error_list",
    "Missed_genotypes",
    "Missed_list",
    "Ratio",
]


def read_whitespace_table(path: Path, expected_header: list[str]) -> list[dict[str, str]]:
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    if not lines:
        raise ValueError(f"{path}: empty table")

    header = lines[0].split()
    if header != expected_header:
        raise ValueError(f"{path}: expected header {expected_header}, found {header}")

    rows = []
    for line_number, line in enumerate(lines[1:], start=2):
        values = line.split()
        if len(values) != len(header):
            raise ValueError(
                f"{path}:{line_number}: expected {len(header)} fields, found {len(values)}"
            )
        rows.append(dict(zip(header, values)))
    return rows


def read_wide_values(path: Path) -> tuple[str, dict[str, str]]:
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    if len(lines) != 2:
        raise ValueError(f"{path}: expected one header and one data row, found {len(lines)} lines")

    header = lines[0].split()
    values = lines[1].split()
    if not header or header[0] != "SAMPLE":
        raise ValueError(f"{path}: first column must be SAMPLE")
    if len(values) != len(header):
        raise ValueError(f"{path}: expected {len(header)} values, found {len(values)}")
    return values[0], dict(zip(header[1:], values[1:]))


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def collate_kir(donor: str, ncopy_dir: Path, genotype_dir: Path, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    copy_rows = read_whitespace_table(
        ncopy_dir / "copy_numbers.txt", ["SAMPLE", "GENE", "COPY_NUMBER"]
    )
    presence_sample, presence = read_wide_values(ncopy_dir / "presence.table.txt")
    depth_sample, depth = read_wide_values(ncopy_dir / "depth_values.txt")
    ratio_sample, ratio = read_wide_values(ncopy_dir / "ratio_values.txt")

    tool_samples = {row["SAMPLE"] for row in copy_rows}
    tool_samples.update([presence_sample, depth_sample, ratio_sample])
    if len(tool_samples) != 1:
        raise ValueError(f"{ncopy_dir}: inconsistent tool sample identifiers: {sorted(tool_samples)}")

    genes = [row["GENE"] for row in copy_rows]
    for name, values in [("presence", presence), ("depth", depth), ("ratio", ratio)]:
        missing = sorted(set(genes) - values.keys())
        extra = sorted(values.keys() - set(genes))
        if missing or extra:
            raise ValueError(f"{ncopy_dir}: {name} genes differ; missing={missing}, extra={extra}")

    copy_number_rows = [
        {
            "donor": donor,
            "gene": row["GENE"],
            "copy_number": row["COPY_NUMBER"],
            "presence": presence[row["GENE"]],
            "depth": depth[row["GENE"]],
            "ratio": ratio[row["GENE"]],
        }
        for row in copy_rows
    ]

    call_paths = sorted(genotype_dir.glob("cds/*/calls/*.calls.txt"))
    report_paths = sorted(genotype_dir.glob("cds/*/reports/*.txt"))
    report_genes = {path.parent.parent.name for path in report_paths}

    call_rows = []
    for path in call_paths:
        gene = path.parent.parent.name
        for row in read_whitespace_table(path, CALL_COLUMNS):
            call_rows.append(
                {
                    "donor": donor,
                    "gene": gene,
                    "copy_number": row["Copy_number"],
                    "calls": row["Calls"],
                    "ratio": row["Ratio"],
                    "missings": row["Missings"],
                    "candidate_report_available": str(gene in report_genes).lower(),
                }
            )

    candidate_rows = []
    for path in report_paths:
        gene = path.parent.parent.name
        for rank, row in enumerate(read_whitespace_table(path, REPORT_COLUMNS), start=1):
            candidate_rows.append(
                {
                    "donor": donor,
                    "gene": gene,
                    "candidate_rank": str(rank),
                    "copy_number": row["Copy_number"],
                    "chromosome": row["Chr"],
                    "allele_a": row["Allele_A"],
                    "allele_b": row["Allele_B"],
                    "allele_c": row["Allele_C"],
                    "allele_d": row["Allele_D"],
                    "tested_genotypes": row["Tested_genotypes"],
                    "valid_genotypes": row["Valid_genotypes"],
                    "error_list": row["Error_list"],
                    "missed_genotypes": row["Missed_genotypes"],
                    "missed_list": row["Missed_list"],
                    "ratio": row["Ratio"],
                }
            )

    copy_number_path = output_dir / f"{donor}.kir.copy_number.tsv"
    calls_path = output_dir / f"{donor}.kir.calls.tsv"
    candidates_path = output_dir / f"{donor}.kir.genotype_candidates.tsv.gz"
    archive_path = output_dir / f"{donor}.kir-mapper.raw.tar.gz"

    write_tsv(
        copy_number_path,
        ["donor", "gene", "copy_number", "presence", "depth", "ratio"],
        copy_number_rows,
    )
    write_tsv(
        calls_path,
        [
            "donor",
            "gene",
            "copy_number",
            "calls",
            "ratio",
            "missings",
            "candidate_report_available",
        ],
        call_rows,
    )
    write_tsv(
        candidates_path,
        [
            "donor",
            "gene",
            "candidate_rank",
            "copy_number",
            "chromosome",
            "allele_a",
            "allele_b",
            "allele_c",
            "allele_d",
            "tested_genotypes",
            "valid_genotypes",
            "error_list",
            "missed_genotypes",
            "missed_list",
            "ratio",
        ],
        candidate_rows,
    )

    # Nextflow stages directory inputs as symlinks by default. Dereference them
    # so the archive contains native kir-mapper files rather than work-dir links.
    with tarfile.open(archive_path, "w:gz", dereference=True) as archive:
        archive.add(ncopy_dir, arcname="ncopy")
        archive.add(genotype_dir, arcname="genotype")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--donor", required=True)
    parser.add_argument("--ncopy-dir", required=True, type=Path)
    parser.add_argument("--genotype-dir", required=True, type=Path)
    parser.add_argument("--output-dir", default=Path("."), type=Path)
    args = parser.parse_args()
    collate_kir(args.donor, args.ncopy_dir, args.genotype_dir, args.output_dir)


if __name__ == "__main__":
    main()

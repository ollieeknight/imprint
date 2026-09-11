#!/usr/bin/env python3
"""Trim and allowlist-correct xGen UDSeq UMIs for DupCaller."""

import argparse
import gzip
import json
from itertools import zip_longest
from pathlib import Path


def open_text(path, mode):
    if not str(path).endswith(".gz"):
        return open(path, mode)
    # Level 9 on the write side dominates this script's runtime (~3-5 MB/s per stream)
    # and these FASTQs are transient work-dir intermediates, not published output.
    return gzip.open(path, mode + "t", **({"compresslevel": 1} if "w" in mode else {}))


def read_fastq(handle):
    while True:
        record = [handle.readline() for _ in range(4)]
        if not record[0]:
            return
        if any(not line for line in record) or not record[0].startswith("@") or not record[2].startswith("+"):
            raise ValueError("Malformed FASTQ record")
        if len(record[1].rstrip("\r\n")) != len(record[3].rstrip("\r\n")):
            raise ValueError("FASTQ sequence and quality lengths differ")
        yield record


def correction(umi, allowlist):
    if len(umi) != 8 or any(base not in "ACGT" for base in umi):
        return None, "invalid"
    if umi in allowlist:
        return umi, "exact"
    matches = [candidate for candidate in allowlist if sum(a != b for a, b in zip(umi, candidate)) == 1]
    return (matches[0], "corrected") if len(matches) == 1 else (None, "unmatched")


def mate_key(name):
    token = name.rstrip().split()[0]
    return token.removesuffix("/1").removesuffix("/2")


def write_record(handle, record, corrected_pair):
    name, sequence, plus, quality = record
    token = name.rstrip().split()[0]
    umi1, umi2 = corrected_pair
    handle.write(f"{token}_{umi1}+{umi2} DB:Z:{umi1}-{umi2}\n")
    handle.write(sequence[8:])
    handle.write(plus)
    handle.write(quality[8:])


def run(args):
    allowlist = {line.strip().upper() for line in Path(args.allowlist).read_text().splitlines() if line.strip()}
    if len(allowlist) != 32 or any(len(umi) != 8 or any(base not in "ACGT" for base in umi) for umi in allowlist):
        raise ValueError("Allowlist must contain exactly 32 unique 8 bp A/C/G/T UMIs")

    counts = {
        "input_pairs": 0,
        "retained_pairs": 0,
        "exact_pairs": 0,
        "corrected_r1": 0,
        "corrected_r2": 0,
        "corrected_both": 0,
        "rejected_invalid": 0,
        "rejected_unmatched": 0,
        "rejected_short": 0,
    }
    with open_text(args.read1, "r") as r1, open_text(args.read2, "r") as r2, \
         open_text(args.output1, "w") as out1, open_text(args.output2, "w") as out2:
        for pair_number, pair in enumerate(zip_longest(read_fastq(r1), read_fastq(r2)), 1):
            rec1, rec2 = pair
            if rec1 is None or rec2 is None:
                raise ValueError("Paired FASTQs contain different record counts")
            if mate_key(rec1[0]) != mate_key(rec2[0]):
                raise ValueError(f"Mate-name mismatch at pair {pair_number}: {rec1[0].strip()} / {rec2[0].strip()}")
            counts["input_pairs"] += 1
            if len(rec1[1].rstrip()) < 8 or len(rec2[1].rstrip()) < 8:
                counts["rejected_short"] += 1
                continue
            umi1, state1 = correction(rec1[1][:8].upper(), allowlist)
            umi2, state2 = correction(rec2[1][:8].upper(), allowlist)
            if umi1 is None or umi2 is None:
                key = "rejected_invalid" if "invalid" in (state1, state2) else "rejected_unmatched"
                counts[key] += 1
                continue
            counts["retained_pairs"] += 1
            if state1 == state2 == "exact":
                counts["exact_pairs"] += 1
            if state1 == "corrected":
                counts["corrected_r1"] += 1
            if state2 == "corrected":
                counts["corrected_r2"] += 1
            if state1 == state2 == "corrected":
                counts["corrected_both"] += 1
            write_record(out1, rec1, (umi1, umi2))
            write_record(out2, rec2, (umi1, umi2))

    Path(args.metrics).write_text(json.dumps(counts, indent=2, sort_keys=True) + "\n")


def self_test():
    allowlist = {"AAAAAAAA", "CCCCCCCC"}
    assert correction("AAAAAAAA", allowlist) == ("AAAAAAAA", "exact")
    assert correction("AAAAAAAT", allowlist) == ("AAAAAAAA", "corrected")
    assert correction("NNNNNNNN", allowlist) == (None, "invalid")
    assert correction("ACGTACGT", allowlist) == (None, "unmatched")
    assert mate_key("@read/1 1:N:0:1\n") == mate_key("@read/2 2:N:0:1\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--read1")
    parser.add_argument("--read2")
    parser.add_argument("--output1")
    parser.add_argument("--output2")
    parser.add_argument("--allowlist")
    parser.add_argument("--metrics")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    missing = [name for name in ("read1", "read2", "output1", "output2", "allowlist", "metrics") if not getattr(args, name)]
    if missing:
        parser.error("missing arguments: " + ", ".join(missing))
    run(args)


if __name__ == "__main__":
    main()

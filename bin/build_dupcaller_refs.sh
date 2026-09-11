#!/usr/bin/env bash
# Build DupCaller repeat, h5 and gene references. Existing outputs are retained.
set -euo pipefail

REF_ROOT=${REF_ROOT:-/sc-projects/sc-proj-cc12-ag-romagnani/ref/imprint}
CONTAINER=${CONTAINER:-/sc-scratch/sc-scratch-cc12-ag-romagnani/apptainer_cache/dupcaller_1.1.2-dev53eb785.sif}
# PERF is a separate pip package; it is not in the DupCaller image.
PERF_PYTHON=${PERF_PYTHON:-python3}
THREADS=${THREADS:-${SLURM_CPUS_PER_TASK:-16}}

FASTA=$REF_ROOT/genome/GRCh38_no_alt_analysis_set.fasta
GTF=$REF_ROOT/characterisation/cnvkit/gencode.v50.basic.annotation.gtf.gz
TARGETS=$REF_ROOT/bed_files/xGen_exome_hyb_panel_v2_hg38/xgen-exome-hyb-panel-v2-targets-hg38.sorted.bed
GERMLINE=$REF_ROOT/gnomAD/gnomad.germline_resource.v4.1.1.hg38.vcf.gz
REPEAT_TSV=$REF_ROOT/udseq/str/GRCh38_perf_repeats.tsv
GENE_BED=$REF_ROOT/udseq/genes/xgen_v2_coding_exons.bed.gz
NOISE_DIR=$REF_ROOT/udseq/noise

log() { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*" >&2; }
die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }
note() { printf '\n=== %s  %s\n' "$1" "$(date -Is)" >>"$REF_ROOT/udseq/PROVENANCE.txt"
         shift; printf '  %s\n' "$@" >>"$REF_ROOT/udseq/PROVENANCE.txt"; }

run_apptainer() {
    local bind_opts=()
    if [[ -n ${APPTAINER_BIND:-} ]]; then
        local clean_bind="${APPTAINER_BIND#--bind }"
        export APPTAINER_BIND="$clean_bind"
        bind_opts=(--bind "$clean_bind")
    elif [[ -d "$REF_ROOT" ]]; then
        bind_opts=(--bind "$REF_ROOT")
    fi
    apptainer exec "${bind_opts[@]}" "$CONTAINER" "$@"
}

# --- repeats -----------------------------------------------------------------
#
# PERF uses `-u 2` to retain short, low-copy repeats. Homopolymers come from FASTA.

do_repeats() {
    [[ -s $REPEAT_TSV ]] && { log "have $REPEAT_TSV"; return; }
    mkdir -p "$(dirname "$REPEAT_TSV")" "$REF_ROOT/udseq"

    "$PERF_PYTHON" -c 'import PERF' 2>/dev/null \
        || die "PERF not importable by $PERF_PYTHON. Install it first: pip install perf_ssr"

    local n_threads="${THREADS:-${SLURM_CPUS_PER_TASK:-4}}"
    log "running PERF: parallel across contigs with $n_threads workers"

    local stage
    stage="$(dirname "$REPEAT_TSV")/.perf_staging"
    rm -rf "$stage"
    mkdir -p "$stage"
    trap 'rm -rf "$stage"' RETURN

    "$PERF_PYTHON" - "$FASTA" "$REPEAT_TSV" "$stage" "$n_threads" <<'PY'
import sys, os, shutil, subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed

fasta_path = sys.argv[1]
output_tsv = sys.argv[2]
stage_dir = sys.argv[3]
threads = max(1, int(sys.argv[4]))

# Split FASTA by contig.
contigs = []
current_ctg = None
current_fh = None

with open(fasta_path, "r") as f:
    for line in f:
        if line.startswith(">"):
            if current_fh:
                current_fh.close()
            current_ctg = line[1:].split()[0]
            contigs.append(current_ctg)
            ctg_fa = os.path.join(stage_dir, f"{current_ctg}.fa")
            current_fh = open(ctg_fa, "w")
            current_fh.write(line)
        else:
            if current_fh:
                current_fh.write(line)
if current_fh:
    current_fh.close()

# Run PERF on one contig.
def run_perf_contig(ctg):
    ctg_fa = os.path.join(stage_dir, f"{ctg}.fa")
    ctg_tsv = os.path.join(stage_dir, f"{ctg}.tsv")
    cmd = [
        sys.executable, "-m", "PERF.core",
        "-m", "1", "-M", "10", "-u", "2",
        "-i", ctg_fa, "-o", ctg_tsv
    ]
    ctg_log = os.path.join(stage_dir, f"{ctg}.log")
    with open(ctg_log, "w") as log_fh:
        res = subprocess.run(cmd, stdout=log_fh, stderr=subprocess.STDOUT)
    if res.returncode != 0:
        with open(ctg_log) as log_fh:
            tail = "".join(log_fh.readlines()[-20:])
        raise RuntimeError(f"PERF failed on {ctg} (exit {res.returncode}):\n{tail}")
    return ctg

# Run the longest contigs first.
sorted_contigs = sorted(
    contigs,
    key=lambda c: os.path.getsize(os.path.join(stage_dir, f"{c}.fa")),
    reverse=True
)

with ThreadPoolExecutor(max_workers=threads) as pool:
    futures = {pool.submit(run_perf_contig, ctg): ctg for ctg in sorted_contigs}
    for fut in as_completed(futures):
        ctg = fut.result()
        fa_to_rm = os.path.join(stage_dir, f"{ctg}.fa")
        if os.path.exists(fa_to_rm):
            os.remove(fa_to_rm)

# Merge PERF rows in FASTA contig order for DupCaller indexing.
tmp_merged = output_tsv + ".tmp"
with open(tmp_merged, "wb") as out_f:
    for ctg in contigs:
        ctg_tsv = os.path.join(stage_dir, f"{ctg}.tsv")
        if not os.path.exists(ctg_tsv):
            continue
        with open(ctg_tsv, "rb") as in_f:
            shutil.copyfileobj(in_f, out_f)

os.replace(tmp_merged, output_tsv)
PY

    rm -rf "$stage"
    log "wrote $REPEAT_TSV ($(wc -l <"$REPEAT_TSV") repeat rows)"
    note repeats "source: $FASTA" "command: PERF.core -m 1 -M 10 -u 2 (parallel $n_threads workers)" "output: $REPEAT_TSV"
}

# --- index -------------------------------------------------------------------

do_index() {
    [[ -s $FASTA.dbs.h5 ]] && { log "have $FASTA.dbs.h5"; return; }
    [[ -s $REPEAT_TSV ]] || die "run the repeats step first"

    # DupCaller writes the h5 beside whichever FASTA it is given. Build through a
    # staging symlink so a failed run cannot destroy a working index.
    local stage=$REF_ROOT/udseq/.index_staging fa
    rm -rf "$stage"; mkdir -p "$stage"
    fa=$stage/$(basename "$FASTA")
    ln -s "$FASTA" "$fa"
    [[ -f $FASTA.fai ]] && ln -s "$FASTA.fai" "$fa.fai"

    log "indexing: single-threaded, hours, ~64GB peak"
    run_apptainer DupCaller.py index -f "$fa" -rt "$REPEAT_TSV"
    check_h5 "$fa" || die "new index failed its checks, left in $stage"

    mv -f "$fa".{ref,tn,hp,str,dbs}.h5 "$(dirname "$FASTA")/"
    rm -rf "$stage"
    note index "reference: $FASTA" "repeat tsv: $REPEAT_TSV" "container: $CONTAINER"
}

# --- genes -------------------------------------------------------------------
#
# DupCaller cuts column 4 at the first underscore to get the gene name
# (Estimate.py: gene_exon.split("_")[0]), so symbols containing one are rewritten.
# Per-gene depth is sum(coverage)/sum(exon lengths), so overlapping exons from
# different transcripts are merged per gene. Exons are clipped to the targets:
# off-panel bases add zero coverage but would still inflate the denominator.

exons_from_gtf() {  # GTF on stdin -> BED4; GTF is 1-based inclusive
    awk -F'\t' 'BEGIN{OFS="\t"}
        $3=="exon" && $9~/gene_type "protein_coding"/ {
            match($9, /gene_name "[^"]+"/); if (!RSTART) next
            gene = substr($9, RSTART+11, RLENGTH-12); gsub(/_/, "-", gene)
            print $1, $4-1, $5, gene
        }'
}

merge_per_gene() {  # input sorted by chrom, gene, start
    awk -F'\t' 'BEGIN{OFS="\t"}
        NR==1 { c=$1; s=$2; e=$3; g=$4; next }
        $1==c && $4==g && $2<=e { if ($3>e) e=$3; next }
        { print c,s,e,g; c=$1; s=$2; e=$3; g=$4 }
        END { if (NR) print c,s,e,g }'
}

number_exons() {
    awk -F'\t' 'BEGIN{OFS="\t"} $4!=g { g=$4; n=0 } { n++; print $1,$2,$3,$4"_"n }'
}

do_genes() {
    [[ -s $GENE_BED ]] && { log "have $GENE_BED"; return; }
    mkdir -p "$(dirname "$GENE_BED")" "$REF_ROOT/udseq"

    gzip -dc "$GTF" | exons_from_gtf | sort -k1,1 -k2,2n \
        | bedtools intersect -a - -b <(sort -k1,1 -k2,2n "$TARGETS") \
        | sort -k1,1 -k4,4 -k2,2n | merge_per_gene \
        | sort -k4,4 -k1,1 -k2,2n | number_exons \
        | sort -k1,1 -k2,2n | bgzip >"$GENE_BED"
    tabix -f -p bed "$GENE_BED"

    log "wrote $GENE_BED ($(gzip -dc "$GENE_BED" | wc -l) exon intervals)"
    note genes "annotation: $GTF" "targets: $TARGETS" \
        "rules: protein_coding exons, clipped to targets, merged per gene" "output: $GENE_BED"
}

# --- verify ------------------------------------------------------------------
# Check h5 contents, BGZF and INFO/AF.

check_h5() {
    run_apptainer python - "$1" <<'PY'
import h5py, sys
fa = sys.argv[1]
want = [f"chr{c}" for c in list(range(1, 23))] + ["chrX"]
fai = {l.split()[0]: int(l.split()[1]) for l in open(fa + ".fai")}
probe, rc = "chr20", 0
# ref/tn/dbs are one value per base; hp is (run length, cut); str is
# (unit length, repeat count, cut).
for suf, ndim in [("ref", 1), ("tn", 1), ("hp", 2), ("str", 2), ("dbs", 1)]:
    with h5py.File(f"{fa}.{suf}.h5", "r") as h:
        d = h[probe]
        missing = [c for c in want if c not in h]
        rows = {"hp": 2, "str": 3}.get(suf)
        ok = (d.ndim == ndim and d.shape[-1] == fai[probe] and not missing
              and (rows is None or d.shape[0] == rows))
        print(f"  {suf}.h5 ndim={d.ndim}/{ndim} len={d.shape[-1]}/{fai[probe]}"
              f" missing={missing or 'none'} {'OK' if ok else 'BAD'}")
        rc |= not ok
        if suf == "str":
            # row 0 is the repeat unit length, painted only from the PERF tsv;
            # all-zero means the index was built without -rt
            n = int((d[0, :] > 0).sum())
            print(f"    STR rows: {n:,} bases binned {'OK' if n else 'BAD (indexed without -rt)'}")
            rc |= not n
sys.exit(rc)
PY
}

do_verify() {
    local rc=0 f noise_files=()
    for prefix in NOISE SNP; do
        if [[ -s $NOISE_DIR/${prefix}_GRCh38.rens.bed.gz ]]; then
            noise_files+=("$NOISE_DIR/${prefix}_GRCh38.rens.bed.gz")
        else
            noise_files+=("$NOISE_DIR/${prefix}.sorted.GRCh38.bed.gz")
        fi
    done
    for f in "$GENE_BED" "${noise_files[@]}"; do
        [[ -s $f && -s $f.tbi ]] || { printf '  MISS  %s (+.tbi)\n' "$f"; rc=1; continue; }
        # pysam fetch needs BGZF; plain gzip passes -s and fails at read time
        file "$f" | grep -q BGZF && printf '  OK    %s\n' "$(basename "$f")" \
            || { printf '  BAD   not BGZF: %s\n' "$f"; rc=1; }
    done

    # call.py swallows a missing AF and treats every alt as AF=1, masking sites
    bcftools view -h "$GERMLINE" | grep -q 'ID=AF,' \
        && printf '  OK    germline INFO/AF\n' \
        || { printf '  BAD   germline VCF has no INFO/AF\n'; rc=1; }

    check_h5 "$FASTA" || rc=1
    [[ $rc -eq 0 ]] && log "all checks passed" || log "checks FAILED"
    return $rc
}

# --- self-test ---------------------------------------------------------------

self_test() {
    # overlapping exons of one gene, a symbol that split("_") would truncate,
    # and a non-coding gene that must be dropped
    local got want
    got=$(printf '%s\n' \
        'chr1	HAVANA	exon	101	200	.	+	.	gene_type "protein_coding"; gene_name "TP53";' \
        'chr1	HAVANA	exon	151	250	.	+	.	gene_type "protein_coding"; gene_name "TP53";' \
        'chr1	HAVANA	exon	401	500	.	+	.	gene_type "protein_coding"; gene_name "TP53";' \
        'chr1	HAVANA	exon	601	700	.	+	.	gene_type "protein_coding"; gene_name "BAD_NAME";' \
        'chr1	HAVANA	exon	801	900	.	+	.	gene_type "lncRNA"; gene_name "NOPE";' \
        | exons_from_gtf | sort -k1,1 -k4,4 -k2,2n | merge_per_gene \
        | sort -k4,4 -k1,1 -k2,2n | number_exons | sort -k1,1 -k2,2n)
    want=$(printf '%s\n' 'chr1	100	250	TP53_1' 'chr1	400	500	TP53_2' \
        'chr1	600	700	BAD-NAME_1' | sort -k1,1 -k2,2n)
    [[ $got == "$want" ]] || { diff <(echo "$want") <(echo "$got"); die "self-test failed (genes)"; }

    echo "self-test OK"
}

# --- main --------------------------------------------------------------------

case ${1:-} in
    --self-test) self_test; exit ;;
    -h|--help) die "usage: $(basename "$0") [--self-test] [str|index|genes|verify ...]" ;;
esac

# unquoted on purpose: quoting collapses the default into a single word
for step in ${@:-repeats index genes verify}; do
    log "=== $step"
    "do_$step"
done

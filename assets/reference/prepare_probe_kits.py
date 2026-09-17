#!/usr/bin/env python3
import os
import argparse

GRCH38_CANONICAL_SIZES = {
    'chr1': 248956422,
    'chr2': 242193529,
    'chr3': 198295559,
    'chr4': 190214555,
    'chr5': 181538259,
    'chr6': 170805979,
    'chr7': 159345973,
    'chr8': 145138636,
    'chr9': 138394717,
    'chr10': 133797422,
    'chr11': 135086622,
    'chr12': 133275309,
    'chr13': 114364328,
    'chr14': 107043718,
    'chr15': 101991189,
    'chr16': 90338345,
    'chr17': 83257441,
    'chr18': 80373285,
    'chr19': 58617616,
    'chr20': 64444167,
    'chr21': 46709983,
    'chr22': 50818468,
    'chrX': 156040895,
    'chrY': 57227415,
    'chrM': 16569
}

WES_CANONICAL_CHRS = {f'chr{i}' for i in range(1, 23)} | {'chrX', 'chrY'}
MITO_CANONICAL_CHRS = {'chrM'}
CHROM_ORDER = {
    **{f'chr{i}': i for i in range(1, 23)},
    'chrX': 23,
    'chrY': 24,
    'chrM': 25,
}

def parse_fai(fai_path):
    if not fai_path:
        print("[*] No FAI supplied. Using built-in GRCh38 sizes.")
        return GRCH38_CANONICAL_SIZES.copy()
    if not os.path.isfile(fai_path):
        raise FileNotFoundError(f"FAI file not found: {fai_path}")
    
    sizes = {}
    print(f"[*] Reading chromosome sizes from: {fai_path}")
    with open(fai_path, encoding='utf-8') as handle:
        for line_number, line in enumerate(handle, start=1):
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 2:
                raise ValueError(f"{fai_path}:{line_number}: expected at least two tab-separated fields")
            try:
                size = int(parts[1])
            except ValueError as error:
                raise ValueError(f"{fai_path}:{line_number}: invalid sequence length '{parts[1]}'") from error
            if size <= 0:
                raise ValueError(f"{fai_path}:{line_number}: sequence length must be positive")
            sizes[parts[0]] = size

    missing = (WES_CANONICAL_CHRS | MITO_CANONICAL_CHRS) - sizes.keys()
    if missing:
        print(f"[*] Warning: {fai_path} is missing canonical contigs: {', '.join(sorted(missing))}")
    return sizes

def read_bed_regions(bed_path, canonical_chrs):
    regions = []
    if not os.path.isfile(bed_path):
        raise FileNotFoundError(f"BED file not found: {bed_path}")

    print(f"[*] Parsing: {bed_path}")
    with open(bed_path, encoding='utf-8') as handle:
        for line_number, line in enumerate(handle, start=1):
            line_str = line.strip()
            if not line_str or line_str.startswith('track') or line_str.startswith('browser') or line_str.startswith('#'):
                continue
            
            parts = line_str.split('\t')
            if len(parts) < 3:
                parts = line_str.split()
                if len(parts) < 3:
                    raise ValueError(f"{bed_path}:{line_number}: expected at least three fields")
            
            chrom = parts[0]
            if chrom in canonical_chrs:
                try:
                    start = int(parts[1])
                    end = int(parts[2])
                except ValueError as error:
                    raise ValueError(f"{bed_path}:{line_number}: invalid BED coordinates") from error
                if start < 0 or end <= start:
                    raise ValueError(f"{bed_path}:{line_number}: expected 0 <= start < end")
                regions.append((chrom, start, end))
    return regions

def pad_and_merge_regions(regions, pad_size, chrom_sizes):
    if not regions:
        return []

    padded = []
    for chrom, start, end in regions:
        if chrom not in chrom_sizes:
            raise ValueError(f"Chromosome size is unavailable for {chrom}")
        max_size = chrom_sizes[chrom]
        if end > max_size:
            raise ValueError(f"Interval {chrom}:{start}-{end} exceeds chromosome length {max_size}")
        new_start = max(0, start - pad_size)
        new_end = min(max_size, end + pad_size)
        padded.append((chrom, new_start, new_end))

    # Group by chromosome to merge
    chrom_regions = {}
    for chrom, start, end in padded:
        if chrom not in chrom_regions:
            chrom_regions[chrom] = []
        chrom_regions[chrom].append((start, end))

    merged_regions = []
    for chrom, intervals in chrom_regions.items():
        intervals.sort(key=lambda x: (x[0], x[1]))
        
        merged = []
        for start, end in intervals:
            if not merged:
                merged.append([start, end])
            else:
                last_start, last_end = merged[-1]
                if start <= last_end:
                    merged[-1][1] = max(last_end, end)
                else:
                    merged.append([start, end])
        
        for start, end in merged:
            merged_regions.append((chrom, start, end))

    return merged_regions

def sort_regions(regions):
    return sorted(regions, key=lambda x: (CHROM_ORDER.get(x[0], 99), x[1], x[2]))

def write_bed(regions, output_path):
    out_dirname = os.path.dirname(output_path)
    if out_dirname:
        os.makedirs(out_dirname, exist_ok=True)
    print(f"[*] Writing to: {output_path}")

    with open(output_path, 'w') as f:
        for chrom, start, end in regions:
            f.write(f"{chrom}\t{start}\t{end}\n")
    print(f"Wrote {len(regions)} regions.")

def process_kit(kit_name, tasks, bed_dir, out_dir, chrom_sizes):
    print(f"\nProcessing {kit_name}")
    for task in tasks:
        input_file = os.path.join(bed_dir, task['input'])
        output_file = os.path.join(out_dir, task['output'])
        
        canonical_chrs = MITO_CANONICAL_CHRS if task.get('is_mito', False) else WES_CANONICAL_CHRS
        
        regions = read_bed_regions(input_file, canonical_chrs)
        
        if not regions:
            raise ValueError(f"No canonical regions found in {input_file}")
        for chrom, start, end in regions:
            if chrom not in chrom_sizes:
                raise ValueError(f"Chromosome size is unavailable for {chrom}")
            if end > chrom_sizes[chrom]:
                raise ValueError(
                    f"Interval {chrom}:{start}-{end} exceeds chromosome length {chrom_sizes[chrom]}"
                )
            
        if task.get('pad_size', 0) > 0:
            print(f"[*] Padding regions by {task['pad_size']}bp and merging.")
            regions = pad_and_merge_regions(regions, task['pad_size'], chrom_sizes)
        
        sorted_regions = sort_regions(regions)
        write_bed(sorted_regions, output_file)

def main():
    parser = argparse.ArgumentParser(description="Prep BED files for imprint pipeline (WES/Mito kits)")
    parser.add_argument("--bed-dir", default=".", help="Base directory containing the 'original' subfolder with raw kits (default: current directory)")
    parser.add_argument("--fai", help="Path to genome fasta .fai index (optional)")
    parser.add_argument("--out-dir", help="Base directory for processed outputs (default: same as --bed-dir)")
    
    args = parser.parse_args()
    
    bed_dir = os.path.abspath(args.bed_dir)
    out_dir = os.path.abspath(args.out_dir) if args.out_dir else bed_dir
    
    original_dir = os.path.join(bed_dir, "original")
    if not os.path.exists(original_dir):
        original_dir = bed_dir
        print(f"[*] 'original' subdirectory not found in {bed_dir}. Reading from {bed_dir} directly.")
    else:
        print(f"[*] Found 'original' subdirectory. Reading raw files from: {original_dir}")
        
    print(f"[*] Processed outputs will be written to: {out_dir}")

    chrom_sizes = parse_fai(args.fai)

    v6_tasks = [
        {
            'input': 'S07604514/S07604514_Regions.bed',
            'output': 'S07604514_hg38/S07604514_Regions.sorted.bed',
            'pad_size': 0
        },
        {
            'input': 'S07604514/S07604514_Regions.bed',
            'output': 'S07604514_hg38/S07604514_Regions.padded50bp.bed',
            'pad_size': 50
        },
        {
            'input': 'S07604514/S07604514_Covered.bed',
            'output': 'S07604514_hg38/S07604514_Covered.bed',
            'pad_size': 0
        }
    ]

    v8_tasks = [
        {
            'input': 'S33266436/S33266436_Regions.bed',
            'output': 'S33266436_hg38/S33266436_Regions.sorted.bed',
            'pad_size': 0
        },
        {
            'input': 'S33266436/S33266436_Regions.bed',
            'output': 'S33266436_hg38/S33266436_Regions.padded50bp.bed',
            'pad_size': 50
        },
        {
            'input': 'S33266436/S33266436_Mergedprobes.bed',
            'output': 'S33266436_hg38/S33266436_MergedProbes.bed',
            'pad_size': 0
        }
    ]

    v7_tasks = [
        {
            'input': 'S31285117/S31285117_Regions.bed',
            'output': 'S31285117_hg38/S31285117_Regions.sorted.bed',
            'pad_size': 0
        },
        {
            'input': 'S31285117/S31285117_Regions.bed',
            'output': 'S31285117_hg38/S31285117_Regions.padded50bp.bed',
            'pad_size': 50
        },
        {
            'input': 'S31285117/S31285117_MergedProbes.bed',
            'output': 'S31285117_hg38/S31285117_MergedProbes.bed',
            'pad_size': 0
        }
    ]

    twist_tasks = [
        {
            'input': 'twist_exome_2.0/hg38_exome_v2.0.2_targets_sorted_validated.re_annotated.bed',
            'output': 'twist_exome_2.0/hg38_exome_v2.0.2_targets_sorted_validated.re_annotated.sorted.bed',
            'pad_size': 0
        },
        {
            'input': 'twist_exome_2.0/hg38_exome_v2.0.2_targets_sorted_validated.re_annotated.bed',
            'output': 'twist_exome_2.0/hg38_exome_v2.0.2_targets_sorted_validated.re_annotated.padded50bp.bed',
            'pad_size': 50
        },
        {
            'input': 'twist_exome_2.0/Twist_MitoPanel_chrM_all_hg38_target.bed',
            'output': 'twist_exome_2.0/Twist_MitoPanel_chrM_all_hg38_target.sorted.bed',
            'pad_size': 0,
            'is_mito': True
        }
    ]

    xgen_exome_v2_tasks = [
        {
            'input': 'xGen_exome_hyb_panel_v2/xgen-exome-hyb-panel-v2-targets-hg38.bed',
            'output': 'xGen_exome_hyb_panel_v2_hg38/xgen-exome-hyb-panel-v2-targets-hg38.sorted.bed',
            'pad_size': 0
        },
        {
            'input': 'xGen_exome_hyb_panel_v2/xgen-exome-hyb-panel-v2-targets-hg38.bed',
            'output': 'xGen_exome_hyb_panel_v2_hg38/xgen-exome-hyb-panel-v2-targets-hg38.padded50bp.bed',
            'pad_size': 50
        },
        {
            'input': 'xGen_exome_hyb_panel_v2/xgen-exome-hyb-panel-v2-probes-hg38.bed',
            'output': 'xGen_exome_hyb_panel_v2_hg38/xgen-exome-hyb-panel-v2-probes-hg38.sorted.bed',
            'pad_size': 0
        }
    ]

    process_kit("Agilent SureSelect v6", v6_tasks, original_dir, out_dir, chrom_sizes)
    process_kit("Agilent SureSelect v7", v7_tasks, original_dir, out_dir, chrom_sizes)
    process_kit("Agilent SureSelect v8", v8_tasks, original_dir, out_dir, chrom_sizes)
    process_kit("Twist Exome 2.0", twist_tasks, original_dir, out_dir, chrom_sizes)
    process_kit("IDT xGen Exome Hyb Panel v2", xgen_exome_v2_tasks, original_dir, out_dir, chrom_sizes)

    print("\nAll probe kits processed.")

if __name__ == "__main__":
    main()

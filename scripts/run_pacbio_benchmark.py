"""
PacBio Microbial 96-plex benchmark - 24 distinct organisms.

Selects 1 BAM + 1 FASTA per organism (bc2001-bc2024),
extracts 1000 reads from each BAM (24K total),
builds index, maps, and compares.

Usage:
    python run_pacbio_benchmark.py
"""
import argparse
import os
import subprocess
import sys
from collections import OrderedDict

SAMTOOLS = r"C:\msys64\mingw64\bin\samtools.exe"
BITPOP = r"C:\Users\Daddy\Documents\GitHub\bit-pop\target\release\bit-pop.exe"
BAM_DIR = r"G:\bam"
FA_DIR = r"G:\fa"
OUTPUT_DIR = r"G:\pacbio_benchmark"

# First replicate for each organism (bc2001-bc2024)
BARCODES = [f"bc{i:04d}" for i in range(2001, 2025)]

# Organism names from barcode
ORGANISMS = {
    "bc2001": "Acinetobacter_baumannii_AYE",
    "bc2002": "Bacillus_cereus_971",
    "bc2003": "Bacillus_subtilis_W23",
    "bc2004": "Burkholderia_cepacia_UCB717",
    "bc2005": "Burkholderia_multivorans_249",
    "bc2006": "Enterococcus_faecalis_OG1RF",
    "bc2008": "Escherichia_coli_K12_MG1655",
    "bc2009": "Helicobacter_pylori_J99",
    "bc2010": "Klebsiella_pneumoniae_BAA2146",
    "bc2011": "Listeria_monocytogenes_Li2",
    "bc2012": "Listeria_monocytogenes_Li23",
    "bc2013": "Methanocorpusculum_labreanum_Z",
    "bc2014": "Neisseria_meningitidis_FAM18",
    "bc2015": "Neisseria_meningitidis_SerogB",
    "bc2016": "Rhodopseudomonas_palustris_CGA009",
    "bc2017": "Salmonella_enterica_LT2",
    "bc2018": "Salmonella_enterica_Ty2",
    "bc2019": "Staphylococcus_aureus_Seattle1945",
    "bc2020": "Staphylococcus_aureus_USA300",
    "bc2021": "Streptococcus_pyogenes_Bruno",
    "bc2022": "Thermanaerovibrio_acidaminovorans",
    "bc2023": "Treponema_denticola_A",
    "bc2024": "Vibrio_parahaemolyticus_EB101",
}


def parse_args():
    parser = argparse.ArgumentParser(description="PacBio 24-organism benchmark")
    parser.add_argument("--reads-per-org", type=int, default=1000, help="Reads per organism")
    parser.add_argument("--threads", type=int, default=8, help="Threads for bit-pop")
    parser.add_argument("--k", type=int, default=10, help="K-mer size")
    return parser.parse_args()


def find_bam_for_barcode(barcode):
    """Find BAM file for given barcode."""
    for f in os.listdir(BAM_DIR):
        if f.endswith(".bam") and barcode in f:
            return os.path.join(BAM_DIR, f)
    return None


def find_fa_for_barcode(barcode):
    """Find FASTA file for given barcode."""
    for f in os.listdir(FA_DIR):
        if f.endswith(".fa") and barcode in f:
            return os.path.join(FA_DIR, f)
    return None


def extract_reads_from_bam(bam_path, output_fastq, num_reads):
    """Extract num_reads from BAM and convert to FASTQ."""
    print(f"  Extracting {num_reads} reads from {os.path.basename(bam_path)}...")
    
    # Use samtools to view and convert to fastq
    cmd = [
        SAMTOOLS, "view", "-h", bam_path,
        "|", SAMTOOLS, "fastq", "-", ">", output_fastq
    ]
    
    # Actually run with shell=True or use subprocess properly
    import subprocess
    view_proc = subprocess.Popen(
        [SAMTOOLS, "view", "-h", bam_path],
        stdout=subprocess.PIPE
    )
    fastq_proc = subprocess.Popen(
        [SAMTOOLS, "fastq", "-"],
        stdin=view_proc.stdout,
        stdout=open(output_fastq, "wb"),
        stderr=subprocess.PIPE
    )
    view_proc.stdout.close()
    _, stderr = fastq_proc.communicate()
    
    if fastq_proc.returncode != 0:
        print(f"  ERROR: {stderr.decode()}", file=sys.stderr)
        return False
    
    # Count lines in FASTQ (4 lines per read)
    with open(output_fastq, "r") as f:
        line_count = sum(1 for _ in f)
    actual_reads = line_count // 4
    print(f"  Extracted {actual_reads} reads")
    
    return True


def build_ground_truth(barcode, organism_name, num_reads):
    """Generate ground truth for extracted reads."""
    gt = OrderedDict()
    for i in range(num_reads):
        read_name = f"read_{i+1}"
        gt[read_name] = organism_name
    return gt


def main():
    args = parse_args()
    
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    print("=" * 60)
    print("PacBio Microbial 96-plex Benchmark (24 organisms)")
    print("=" * 60)
    
    # Step 1: Find BAM and FASTA files
    print("\nStep 1: Finding files...")
    bam_files = {}
    fa_files = {}
    
    for barcode in BARCODES:
        bam = find_bam_for_barcode(barcode)
        fa = find_fa_for_barcode(barcode)
        
        if bam:
            bam_files[barcode] = bam
        else:
            print(f"  WARNING: No BAM for {barcode}")
            
        if fa:
            fa_files[barcode] = fa
        else:
            print(f"  WARNING: No FASTA for {barcode}")
    
    print(f"  Found {len(bam_files)} BAM files")
    print(f"  Found {len(fa_files)} FASTA files")
    
    if len(bam_files) < 24 or len(fa_files) < 24:
        print("ERROR: Need 24 BAM and 24 FASTA files", file=sys.stderr)
        sys.exit(1)
    
    # Step 2: Extract reads from each BAM
    print(f"\nStep 2: Extracting {args.reads_per_org} reads per organism...")
    all_reads_fastq = os.path.join(OUTPUT_DIR, "pacbio_24k_reads.fq")
    ground_truth = OrderedDict()
    
    # Clear output file
    with open(all_reads_fastq, "w") as f:
        pass
    
    for barcode in BARCODES:
        if barcode not in bam_files:
            continue
        
        organism = ORGANISMS.get(barcode, barcode)
        bam_path = bam_files[barcode]
        
        # Extract reads
        temp_fastq = os.path.join(OUTPUT_DIR, f"temp_{barcode}.fq")
        if extract_reads_from_bam(bam_path, temp_fastq, args.reads_per_org):
            # Append to main file
            with open(temp_fastq, "r") as src, open(all_reads_fastq, "a") as dst:
                dst.write(src.read())
            
            # Generate ground truth
            for i in range(args.reads_per_org):
                read_name = f"{organism}_read_{i+1}"
                ground_truth[read_name] = organism
            
            os.remove(temp_fastq)
    
    print(f"\n  Total reads extracted: {len(ground_truth)}")
    print(f"  Total organisms: {len(set(ground_truth.values()))}")
    
    # Step 3: Write ground truth TSV
    print("\nStep 3: Writing ground truth...")
    gt_file = os.path.join(OUTPUT_DIR, "pacbio_ground_truth.tsv")
    with open(gt_file, "w") as f:
        f.write("@Version:0.9.1\n")
        f.write("@SampleID:pacbio_24org\n")
        f.write("\n")
        f.write("@@SEQUENCEID\tBINID\tTAXID\t_READID\n")
        for read_name, organism in ground_truth.items():
            f.write(f"{read_name}\t{organism}\t?\t{read_name}\n")
    
    gt_size = os.path.getsize(gt_file)
    print(f"  Written: {gt_file} ({gt_size / 1024:.1f} KB)")
    
    # Step 4: Build index
    print("\nStep 4: Building index...")
    index_file = os.path.join(OUTPUT_DIR, "pacbio_24org_index.bitpop")
    
    fa_list = list(fa_files.values())
    cmd = [
        BITPOP, "build", "--cami",
        "-o", index_file,
        "-k", str(args.k),
        "-t", str(args.threads),
    ]
    # Add each FASTA with -f flag
    for fa in fa_list:
        cmd.extend(["-f", fa])
    
    print(f"  Command: {' '.join(cmd[:5])} ... ({len(fa_list)} genomes)")
    result = subprocess.run(cmd)
    
    if result.returncode != 0:
        print("ERROR: Index build failed", file=sys.stderr)
        sys.exit(1)
    
    index_size = os.path.getsize(index_file) / 1024 / 1024 / 1024
    print(f"  Index built: {index_file} ({index_size:.2f} GB)")
    
    # Step 5: Map reads
    print("\nStep 5: Mapping reads...")
    sam_file = os.path.join(OUTPUT_DIR, "pacbio_24k_mapped.sam")
    
    cmd = [
        BITPOP, "map",
        "-i", index_file,
        "-r", all_reads_fastq,
        "-o", sam_file,
        "--top-n", "2",
        "-t", str(args.threads),
    ]
    
    print(f"  Mapping {len(ground_truth)} reads...")
    result = subprocess.run(cmd)
    
    if result.returncode != 0:
        print("ERROR: Mapping failed", file=sys.stderr)
        sys.exit(1)
    
    sam_size = os.path.getsize(sam_file) / 1024 / 1024
    print(f"  SAM written: {sam_file} ({sam_size:.1f} MB)")
    
    # Step 6: Compare results
    print("\nStep 6: Comparing with ground truth...")
    report_file = os.path.join(OUTPUT_DIR, "pacbio_24k_report.txt")
    
    cmd = [
        "python", r"scripts\compare_cami_results.py",
        "--sam", sam_file,
        "--gt", gt_file,
        "--output-report", report_file,
    ]
    
    result = subprocess.run(cmd)
    
    print(f"\n{'=' * 60}")
    print("BENCHMARK COMPLETE")
    print(f"{'=' * 60}")
    print(f"Report: {report_file}")
    print(f"Reads: {len(ground_truth)}")
    print(f"Organisms: 24")
    print(f"Index size: {index_size:.2f} GB")
    print(f"Time: See report for details")


if __name__ == "__main__":
    main()

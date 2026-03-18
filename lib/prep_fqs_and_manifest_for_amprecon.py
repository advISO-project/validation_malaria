#!/usr/bin/env python3
import argparse
import os
import gzip

# ----------------------------
# Argument parsing
# ----------------------------
parser = argparse.ArgumentParser(
    description="Generate interleaved FASTQ files from paired gzipped inputs "
                "and produce a manifest CSV."
)
parser.add_argument(
    "--insdc_manifest",
    required=True,
    help="Path to input CSV manifest (with columns: sample,sample_acc,library,experiment_acc,run_acc)"
)
parser.add_argument(
    "--insdc_fastq_folder",
    required=True,
    help="Path to folder containing paired gzipped FASTQ files named <run_acc>_1.fastq.gz and <run_acc>_2.fastq.gz"
)
parser.add_argument(
    "--output_folder",
    required=True,
    help="Folder where interleaved FASTQs and output manifest CSV will be written"
)

args = parser.parse_args()

input_csv = args.insdc_manifest
input_path = args.insdc_fastq_folder
output_folder = args.output_folder

os.makedirs(output_folder, exist_ok=True)
output_csv = os.path.join(output_folder, "amprecon_manifest.csv")

# ----------------------------
# Mapping library suffix → primer panel
# ----------------------------
panel_map = {
    "GRC1": "PFA_GRC1_v1.0",
    "GRC2": "PFA_GRC2_v1.0",
    "SPEC": "PFA_SPEC"
}

# ----------------------------
# Function to interleave paired FASTQs
# ----------------------------
def interleave_fastq(run_acc, input_path, output_folder):
    file1 = os.path.join(input_path, f"{run_acc}_1.fastq.gz")
    file2 = os.path.join(input_path, f"{run_acc}_2.fastq.gz")
    output_file = os.path.join(output_folder, f"{run_acc}.fastq")
    
    if not os.path.exists(file1) or not os.path.exists(file2):
        print(f"Warning: missing pair for {run_acc}, skipping.")
        return None
    
    with gzip.open(file1, "rt") as f1, gzip.open(file2, "rt") as f2, open(output_file, "w") as out:
        while True:
            r1 = [f1.readline() for _ in range(4)]
            r2 = [f2.readline() for _ in range(4)]
            if not r1[0] or not r2[0]:
                break
            out.writelines(r1)
            out.writelines(r2)
    
    return output_file

# ----------------------------
# Process CSV and generate outputs
# ----------------------------
with open(input_csv, "r") as fin, open(output_csv, "w") as fout:
    fout.write("sample_id,primer_panel,fastq_path\n")
    next(fin)  # skip header
    
    for line in fin:
        line = line.strip()
        if not line:
            continue
        sample, sample_acc, library, experiment_acc, run_acc = line.split(",")
        suffix = library.split("_")[-1]
        primer_panel = panel_map.get(suffix, "UNKNOWN")
        
        interleaved_file = interleave_fastq(run_acc, input_path, output_folder)
        if interleaved_file:
            fout.write(f"{sample},{primer_panel},{interleaved_file}\n")

print(f"All done! Interleaved FASTQs and manifest CSV written to {output_folder}")
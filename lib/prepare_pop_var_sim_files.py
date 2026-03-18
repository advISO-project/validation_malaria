#!/usr/bin/env python3

import os
import csv
import json
import argparse
from collections import defaultdict

from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Seq import Seq


# -------------------------
# Genetic code helpers
# -------------------------

table = CodonTable.unambiguous_dna_by_name["Standard"]

def translate_codon(codon):
    return table.forward_table.get(codon.upper(), "*")

def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())


# -------------------------
# Parse genome & GFF
# -------------------------

def load_genome(fasta):
    return SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))

def load_cds_map(gff_file):
    """
    Simplified GFF3 parser: collect CDS coordinates per gene using the 'gene_id' attribute.
    Returns: dict of gene_id -> (chrom, strand, list of genomic positions in CDS order)
    """
    from collections import defaultdict

    cds_dict = defaultdict(list)
    strand_dict = {}
    chrom_dict = {}

    with open(gff_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            cols = line.split("\t")
            chrom, source, feature_type, start, end, score, strand, phase, attributes = cols
            start, end = int(start), int(end)

            if feature_type != "CDS":
                continue

            # parse attributes into a dict
            attr_dict = {}
            for attr in attributes.split(";"):
                if "=" in attr:
                    key, val = attr.split("=", 1)
                    attr_dict[key] = val

            gene_id = attr_dict.get("gene_id")
            if not gene_id:
                continue  # skip CDS without gene_id

            # record chrom and strand
            strand_dict[gene_id] = strand
            chrom_dict[gene_id] = chrom

            # add all positions for this CDS
            cds_dict[gene_id].extend(range(start, end + 1))

    # final CDS map
    cds_map = {}
    for gene_id, positions in cds_dict.items():
        if strand_dict[gene_id] == "-":
            positions.reverse()
        cds_map[gene_id] = (chrom_dict[gene_id], strand_dict[gene_id], positions)

    return cds_map


# -------------------------
# Codon → SNP proposal
# -------------------------

def propose_alt_codon(ref_codon, target_aa):
    if translate_codon(ref_codon) == target_aa:
        return ref_codon

    bases = ["A","C","G","T"]

    # try single SNP
    for i in range(3):
        for b in bases:
            if b == ref_codon[i]:
                continue
            new = ref_codon[:i] + b + ref_codon[i+1:]
            if translate_codon(new) == target_aa:
                return new

    raise ValueError(f"Cannot convert {ref_codon} → {target_aa}")


def codon_positions(cds_coords, strand, aa_pos):
    codon_index = aa_pos - 1
    codon = cds_coords[codon_index*3 : codon_index*3 + 3]

    if strand == "-":
        codon = codon[::-1]

    return codon


# -------------------------
# Parse haplotype CSV
# -------------------------

def split_allele(val):
    if "/" in val:
        left, right = val.split("/")
        return left, right
    return val, val


def build_haplotypes(row, markers):
    hap1 = {}
    hap2 = {}

    for m in markers:
        a1, a2 = split_allele(row[m])
        hap1[m] = a1
        hap2[m] = a2

    if hap1 == hap2:
        return [hap1]

    return [hap1, hap2]


# -------------------------
# SNP extraction
# -------------------------

def haplotype_snps(hap, gene_map, genome, cds_map):
    snps = []

    for marker, aa in hap.items():
        gene, pos = marker.split("_")
        pos = int(pos)

        gene_id = gene_map[gene]
        chrom, strand, cds_coords = cds_map[gene_id]

        positions = codon_positions(cds_coords, strand, pos)

        ref_codon = "".join(str(genome[chrom].seq[p-1]) for p in positions)

        if strand == "-":
            ref_codon = reverse_complement(ref_codon)

        ref_aa = translate_codon(ref_codon)

        if ref_aa == aa:
            continue

        alt_codon = propose_alt_codon(ref_codon, aa)

        if strand == "-":
            ref_codon = reverse_complement(ref_codon)
            alt_codon = reverse_complement(alt_codon)

        for i in range(3):
            if ref_codon[i] != alt_codon[i]:
                snps.append((chrom, positions[i], ref_codon[i], alt_codon[i]))

    return snps


def write_vcf_like(path, hap_name, snps):
    with open(path, "w") as f:
        f.write("#CHROM\tPOS\tREF\tALT\n")

        for s in snps:
            f.write("\t".join(map(str, s)) + "\n")


# -------------------------
# Panel manifest
# -------------------------

def write_panel_manifest(vcf_files, bed_files, genome_file, out_csv):
    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["haplotype","base_fasta","variants_file","bed_file"])

        for hap, vcf in vcf_files.items():
            for suffix, bed in bed_files.items():
                writer.writerow([
                    f"{hap}_{suffix}",
                    genome_file,
                    vcf,
                    bed
                ])


# -------------------------
# JSON genotype output
# -------------------------

def write_json(samples, out_json, reads):
    """
    samples: dict mapping sample name -> list of haplotype IDs
    out_json: path to output JSON
    reads: dict with keys GRC1, GRC2, SPEC
    """
    entries = []

    for sample, genotypes in samples.items():
        if len(genotypes) == 1:
            proportions = [1.0]
        else:
            proportions = [0.5, 0.5]

        for panel in ["GRC1", "GRC2", "SPEC"]:
            g_suf = "_" + panel 
            entries.append({
                "sample_id": f"{sample}_{panel}",
                "genotypes": [g + g_suf for g in genotypes],
                "proportions": proportions,
                "num_reads": reads[panel]
            })

    with open(out_json, "w") as f:
        json.dump(entries, f, indent=2)

# -------------------------
# MAIN
# -------------------------

def main():
    p = argparse.ArgumentParser()

    p.add_argument("--fasta", required=True)
    p.add_argument("--gff", required=True)
    p.add_argument("--csv", required=True)

    p.add_argument("--bed_grc1", required=True)
    p.add_argument("--bed_grc2", required=True)
    p.add_argument("--bed_spec", required=True)

    p.add_argument("--reads_grc1", type=int, required=True)
    p.add_argument("--reads_grc2", type=int, required=True)
    p.add_argument("--reads_spec", type=int, required=True)

    args = p.parse_args()

    genome = load_genome(args.fasta)
    cds_map = load_cds_map(args.gff)

    gene_map = {
        "dhfr": "PF3D7_0417200",
        "dhps": "PF3D7_0810800",
        "crt": "PF3D7_0709000",
        "kelch13": "PF3D7_1343700"
    }

    bed_files = {
        "GRC1": args.bed_grc1,
        "GRC2": args.bed_grc2,
        "SPEC": args.bed_spec
    }

    reads = {
        "GRC1": args.reads_grc1,
        "GRC2": args.reads_grc2,
        "SPEC": args.reads_spec
    }

    vcf_files = {}
    sample_genotypes = {}

    with open(args.csv) as f:
        reader = csv.DictReader(f)

        markers = [c for c in reader.fieldnames if "_" in c and c.split("_")[0] in gene_map]

        for row in reader:
            sample = row["sample"]

            haplotypes = build_haplotypes(row, markers)

            genotype_ids = []

            for i, hap in enumerate(haplotypes):
                hap_name = sample if len(haplotypes)==1 else f"{sample}_h{i+1}"

                snps = haplotype_snps(hap, gene_map, genome, cds_map)

                vcf_name = hap_name + ".vcf"
                vcf_path = os.path.join(".", vcf_name)
                write_vcf_like(vcf_path, hap_name, snps)

                vcf_files[hap_name] = vcf_name
                genotype_ids.append(hap_name)

            sample_genotypes[sample] = genotype_ids

    write_panel_manifest(
        vcf_files,
        bed_files,
        args.fasta,
        os.path.join(".", "haplotype_panel_manifest.csv")
    )

    write_json(
        sample_genotypes,
        os.path.join(".", "sample_design.json"),
        reads
    )


if __name__ == "__main__":
    main()

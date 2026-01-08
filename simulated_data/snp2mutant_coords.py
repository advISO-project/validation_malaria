#!/usr/bin/env python3
import sys
import re
from collections import defaultdict

def parse_attributes(attr_str):
    """Parse GFF attribute column into dict."""
    attrs = {}
    for part in attr_str.strip().split(";"):
        if not part:
            continue
        if "=" in part:
            k,v = part.split("=",1)
            attrs[k] = v
        elif " " in part:
            # fallback for space-separated key value
            k,v = part.split(" ",1)
            attrs[k] = v
    return attrs

def find_cds_for_gene(gff_path, query):
    """
    Returns:
      chrom (seqid), strand (+1 or -1), list of (start,end) tuples for CDS in transcript order,
      matched_ids -> list of attribute IDs found (for debugging)
    """
    # Map parent/id -> list of CDS tuples and seqid/strand
    cds_by_parent = defaultdict(list)
    # Track any attributes seen that contain query (for diagnostics)
    matching_attrs = set()
    # Also track lines where gene/transcript types include query
    candidate_lines = []

    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            seqid, source, ftype, start_s, end_s, score, strand_s, phase, attr_s = cols
            start = int(start_s)
            end = int(end_s)
            attrs = parse_attributes(attr_s)

            # Check common attribute keys for matching
            for key in ("ID","Parent","Name","gene_id","transcript_id","Dbxref"):
                val = attrs.get(key)
                if val and query in val:
                    matching_attrs.add(f"{key}={val}")
                    candidate_lines.append((seqid, ftype, start, end, strand_s, attrs))

            # If it's a CDS, record under its Parent (often "Parent=PF3D7_0709000.1")
            if ftype == "CDS":
                parent = attrs.get("Parent") or attrs.get("Gene") or attrs.get("Transcript")
                # Parent can be comma-separated
                if parent:
                    for p in parent.split(","):
                        cds_by_parent[p].append((seqid, start, end, strand_s))
                else:
                    # fallback: sometimes ID is used for CDS
                    cid = attrs.get("ID")
                    if cid:
                        cds_by_parent[cid].append((seqid, start, end, strand_s))

    # If we found no matching_attrs, try scanning names that contain the short query substring
    if not matching_attrs:
        # Look for any parent id that contains the query substring
        candidates = [p for p in cds_by_parent.keys() if query in p]
        if candidates:
            for p in candidates:
                matching_attrs.add(f"Parent_contains={p}")

    # Build final CDS list by selecting parents that match
    chosen = None
    cds_list = []
    chrom = None
    strand = None
    for parent, cds_entries in cds_by_parent.items():
        if query in parent or any(query in m for m in matching_attrs):
            # prefer parent matches
            if query in parent:
                chosen = parent
                cds_list = cds_entries
                break

    # If not chosen, but there are matching_attrs that point to IDs, try to pick one
    if chosen is None:
        # pick any parent that contains the query substring
        for parent in cds_by_parent:
            if query in parent:
                chosen = parent
                cds_list = cds_by_parent[parent]
                break

    # If still no cds_list, try any parent where the ID startswith query (useful for .1 suffix)
    if not cds_list:
        for parent in cds_by_parent:
            if parent.startswith(query):
                chosen = parent
                cds_list = cds_by_parent[parent]
                break

    # If still empty, return diagnostics
    if not cds_list:
        return None, None, None, sorted(matching_attrs)

    # Ensure seqid and strand are consistent
    seqids = set(e[0] for e in cds_list)
    strands = set(e[3] for e in cds_list)
    if len(seqids) > 1 or len(strands) > 1:
        raise ValueError(f"CDS entries for {chosen} have inconsistent seqid/strand: seqids={seqids}, strands={strands}")

    chrom = cds_list[0][0]
    strand = +1 if cds_list[0][3] == "+" else -1

    # convert to (start,end) and sort in transcript order
    exons = [(s,e) for (_,s,e,_) in cds_list]
    if strand == 1:
        exons.sort(key=lambda x: x[0])   # ascending start
    else:
        exons.sort(key=lambda x: x[0], reverse=True)  # reverse genomic order for transcript order

    # Normalize as integer inclusive coords
    exons_norm = [(int(a), int(b)) for a,b in exons]
    return chrom, strand, exons_norm, sorted(matching_attrs)

def aa_to_genomic_from_exons(exons, aa_pos, strand):
    """
    exons: list of (start,end) in transcript order (1-based inclusive)
    aa_pos: 1-based amino acid position
    returns (genomic_start, genomic_end) inclusive (1-based)
    """
    # codon positions in CDS (1-based)
    codon_start_cds = (aa_pos - 1) * 3 + 1
    codon_end_cds = codon_start_cds + 2

    # iterate exons accumulating length
    running = 0
    g_positions = []
    for start,end in exons:
        exon_len = end - start + 1
        exon_cds_start = running + 1
        exon_cds_end = running + exon_len
        # overlap?
        overlap_start = max(exon_cds_start, codon_start_cds)
        overlap_end = min(exon_cds_end, codon_end_cds)
        if overlap_start <= overlap_end:
            # compute genomic coords for that overlap
            offset_start = overlap_start - exon_cds_start  # 0-based offset into exon
            offset_end = overlap_end - exon_cds_start
            if strand == 1:
                g_s = start + offset_start
                g_e = start + offset_end
            else:
                # transcript order for negative strand: exons were listed in reverse genomic order,
                # but start/end are still genomic coordinates where start < end.
                # For negative strand, the first base of the exon in transcript is exon_end (genomic)
                g_s = end - offset_start
                g_e = end - offset_end
                # ensure g_s <= g_e
                g_s, g_e = min(g_s,g_e), max(g_s,g_e)
            g_positions.append((g_s, g_e))
        running += exon_len

    if not g_positions:
        raise ValueError("Codon not found within exons (check aa_pos and exon boundaries).")

    # merge contiguous ranges (should normally be 1 range, except codon spanning exon junction -> 2 ranges)
    g_positions.sort()
    merged_start = g_positions[0][0]
    merged_end = g_positions[0][1]
    for s,e in g_positions[1:]:
        if s <= merged_end + 1:
            merged_end = max(merged_end, e)
        else:
            # non-contiguous (rare for codon)
            pass

    return merged_start, merged_end

def sequence_for_pos(fasta_file, chrom, strand, start, end):
    """
    fasta_file: path (str) to genome FASTA file
    chrom: name of chromosome, must match (case insensitive) a valid ID from FASTA file (first word in the header)
    strand: strand (int, 1:'+', -1:'-')
    start: 1-based start pos
    end: 1-based end pos
    returns NA sequence at given position (str)
    NOTE: this is a very simple FASTA parser that doesn't do any error handling
    """
    seek=False
    seq=''
    with open(fasta_file) as fh:
        for line in fh:
            line=line.strip()
            if line.startswith('>'):
                seek=False
                seq_id=re.match(r'^> *(\S+).*', line)[1]
                if seq_id.lower()==chrom.lower():
                    seek=True
            elif seek:
                seq=seq+line
                if len(seq)>=end:
                    if strand==1:
                        return seq[start-1:end]
                    else:
                        return _revcom(seq[start-1:end])
    raise ValueError('could not find chromosome or requested position in chromosome')   
                
def _revcom(seq):
    """
    seq: sequence string
    returns: reverse complement of sequence (uppercased)
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement.get(base, base) for base in reversed(seq.upper()))

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python aa_to_genomic.py <PlasmoDB-68.gff> <GENE_ID_query> <AA_position>")
        print("Example: python aa_to_genomic.py PlasmoDB-68_Pfalciparum3D7.gff PF3D7_0709000 76")
        sys.exit(1)

    gff = sys.argv[1]
    geneq = sys.argv[2].strip()
    aa = int(sys.argv[3])

    chrom, strand, exons, diagnostics = find_cds_for_gene(gff, geneq)

    if exons is None:
        print("No CDS found for your query. Diagnostics / matching attributes found in GFF:")
        for m in diagnostics:
            print("  ", m)
        # Also give a hint to search file for the geneq substring
        print("\nTip: try running this to find candidate lines in the GFF that contain your query:")
        print(f"  grep -n '{geneq}' {gff} | head -n 50")
        sys.exit(2)

    print(f"Found CDS on {chrom}, strand {'+' if strand==1 else '-'}, {len(exons)} exon(s).")
    print("Exons (transcript order):")
    for s,e in exons:
        print(f"  {s}-{e}  (length {e-s+1})")

    g_s, g_e = aa_to_genomic_from_exons(exons, aa, strand)
    print(f"\nAmino acid {aa} -> codon genomic coordinates (1-based inclusive):")
    print(f"Chromosome: {chrom}")
    print(f"Genomic codon: {chrom}:{g_s}-{g_e}  (strand {'+' if strand==1 else '-'})")

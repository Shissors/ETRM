import subprocess
import argparse
import os
import random
import pandas as pd
import json
from collections import Counter
from Bio import AlignIO, Phylo, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

# -------------------- MAFFT Alignment --------------------
def run_mafft(input_fasta, output_fasta="aligned.fasta", mafft_path="mafft", options="--auto"):
    if not os.path.exists(input_fasta):
        raise FileNotFoundError(f"Input file not found: {input_fasta}")

    cmd = [mafft_path] + options.split() + [input_fasta]
    print(f"\n⚡ Running MAFFT command: {' '.join(cmd)}")
    print("This may take a while for large alignments...")

    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        print("❌ MAFFT failed:")
        print(result.stderr)
        raise RuntimeError("MAFFT alignment failed.")

    with open(output_fasta, "w") as f:
        f.write(result.stdout)

    print(f"✅ Aligned sequences saved to: {output_fasta}")
    return output_fasta

# -------------------- Ancestral Sequence Reconstruction --------------------
def reconstruct_weighted_ancestor(alignment_file, output_file="Consensus_Sequence.fasta",
                                  ancestor_id="Pigeon_HBB_Ancestor", tree_file="calculated_tree.xml",
                                  newick_file="calculated_tree.nwk"):
    print("\n🔬 Starting ancestral sequence reconstruction...")
    alignment = AlignIO.read(alignment_file, "fasta")

    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(alignment)
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)
    tree.root_at_midpoint()

    seq_length = alignment.get_alignment_length()
    ancestral_seq = []

    distance_map = {clade.name: tree.root.distance(clade) for clade in tree.get_terminals()}
    max_dist = max(distance_map.values()) if distance_map else 1.0
    if max_dist == 0:
        max_dist = 1.0

    for i in range(seq_length):
        column_data = Counter()
        for rec in alignment:
            base = rec.seq[i]
            try:
                distance = distance_map.get(rec.id, max_dist)
                normalized_distance = distance / max_dist
                weight = 1.0 / (normalized_distance + 0.01)
            except Exception:
                weight = 1.0
            column_data[base] += weight

        most_likely_base = column_data.most_common(1)[0][0] if column_data else '?'
        ancestral_seq.append(most_likely_base)

    final_seq = "".join(ancestral_seq)
    rec = SeqRecord(Seq(final_seq), id=ancestor_id,
                    description="Reconstruction using Neighbor-Joining tree weights.")
    SeqIO.write([rec], output_file, "fasta")

    print("\n✅ --- Result Summary ---")
    print(f"Reconstruction Method: Neighbor-Joining Tree Weighted Consensus")
    print(f"Reconstructed sequence length: {len(final_seq)}")
    print(f"Sequence written to: {output_file}")

    # Save PhyloXML tree
    try:
        Phylo.write(tree, tree_file, "phyloxml")
        print(f"Tree saved to {tree_file}")
    except Exception as e:
        print(f"Could not save tree file: {e}")

    # Export Newick tree to file
    export_newick_tree(tree_file, newick_file)

def export_newick_tree(tree_file="calculated_tree.xml", newick_file="calculated_tree.nwk"):
    print(f"\n📂 Exporting Newick format from '{tree_file}'...")
    if not os.path.exists(tree_file):
        print(f"Error: Tree file '{tree_file}' not found. Skipping Newick export.")
        return
    try:
        tree = Phylo.read(tree_file, "phyloxml")
        newick_string = tree.format('newick')
        with open(newick_file, "w") as f:
            f.write(newick_string)
        print(f"Newick tree saved to: {newick_file}")
    except Exception as e:
        print(f"An error occurred while processing the tree file: {e}")

# -------------------- Protein Mutation Simulation --------------------
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
AMINO_ACID_CLASSES = {
    'A': 1, 'V': 1, 'I': 1, 'L': 1, 'M': 1, 'F': 1, 'W': 1, 'P': 1,
    'G': 1,
    'S': 2, 'T': 2, 'C': 2, 'N': 2, 'Q': 2, 'Y': 2, 'H': 2,
    'D': 3, 'E': 3, 'K': 3, 'R': 3
}

def is_lethal_change(old_aa, new_aa, position_weight, LETHAL_CHANGE_REJECTION_RATE):
    old_class = AMINO_ACID_CLASSES.get(old_aa, 0)
    new_class = AMINO_ACID_CLASSES.get(new_aa, 0)
    property_change = abs(old_class - new_class)
    if position_weight < 0.2 and property_change >= 2:
        return random.random() < LETHAL_CHANGE_REJECTION_RATE
    if 0.2 <= position_weight < 0.6 and property_change >= 2:
        return random.random() < LETHAL_CHANGE_REJECTION_RATE
    return False

def mutate_residue_weighted(old_aa):
    old_class = AMINO_ACID_CLASSES.get(old_aa, 0)
    weights = []
    candidates = list(AMINO_ACIDS)
    for new_aa in candidates:
        if new_aa == old_aa:
            weights.append(0)
            continue
        new_class = AMINO_ACID_CLASSES.get(new_aa, 0)
        property_change = abs(old_class - new_class)
        if property_change == 0:
            weights.append(5)
        elif property_change == 1:
            weights.append(2)
        else:
            weights.append(1)
    return random.choices(candidates, weights=weights, k=1)[0]

def generate_dummy_conservation_weights(length):
    weights = []
    flexible_length = int(length * 0.2)
    conserved_length = int(length * 0.4)
    for i in range(length):
        if i >= (length - conserved_length) // 2 and i < (length + conserved_length) // 2:
            weights.append(0.1)
        elif i < flexible_length or i >= length - flexible_length:
            weights.append(1.0)
        else:
            weights.append(0.5)
    return weights

def simulate_mutations(input_fasta, snapshot_dir, num_snapshots=20,
                       mutation_rate_per_site=0.01, LETHAL_CHANGE_REJECTION_RATE=0.95,
                       combined_fasta_file="all_snapshots.fasta"):
    if not os.path.exists(snapshot_dir):
        os.makedirs(snapshot_dir)

    try:
        record = SeqIO.read(input_fasta, "fasta")
        sequence = str(record.seq).upper().replace('*', '').replace('-', '')
    except Exception as e:
        print(f"Error reading {input_fasta}: {e}")
        return

    seq_length = len(sequence)
    CONSERVATION_WEIGHTS = generate_dummy_conservation_weights(seq_length)
    mutation_log = []
    current_seq = list(sequence)
    output_prefix = "protein_snap"
    all_snap_records = []

    print(f"\n🧬 Simulating {num_snapshots} snapshots...")

    for snap in range(num_snapshots):
        new_seq = current_seq.copy()
        mutations_in_snap = 0
        rejected_mutations = 0
        for i in range(seq_length):
            weight = CONSERVATION_WEIGHTS[i]
            effective_rate = mutation_rate_per_site * weight
            if random.random() < effective_rate:
                old = current_seq[i]
                new = mutate_residue_weighted(old)
                if is_lethal_change(old, new, weight, LETHAL_CHANGE_REJECTION_RATE):
                    new_seq[i] = old
                    rejected_mutations += 1
                    mutation_log.append({
                        "snapshot": snap, "position": i+1, "old_aa": old, "new_aa": new,
                        "conservation_weight": weight, "status": "REJECTED (Lethal Change)"
                    })
                else:
                    new_seq[i] = new
                    mutations_in_snap += 1
                    mutation_log.append({
                        "snapshot": snap, "position": i+1, "old_aa": old, "new_aa": new,
                        "conservation_weight": weight, "status": "ACCEPTED"
                    })

        snap_seq_str = "".join(new_seq)
        snap_record = SeqRecord(Seq(snap_seq_str),
                                id=f"{output_prefix}_{snap:02d}",
                                description=f"Protein Snapshot {snap}")
        all_snap_records.append(snap_record)
        current_seq = new_seq
        print(f"Snapshot {snap:02d}: {mutations_in_snap} ACCEPTED, {rejected_mutations} REJECTED")

    # Write all snapshots to a single FASTA
    combined_fasta_path = os.path.join(snapshot_dir, combined_fasta_file)
    SeqIO.write(all_snap_records, combined_fasta_path, "fasta")
    print(f"\n✅ All snapshots written to: {combined_fasta_path}")

    # Save mutation log
    csv_filepath = os.path.join(snapshot_dir, "amino_acid_mutation_log.csv")
    pd.DataFrame(mutation_log).to_csv(csv_filepath, index=False)
    print(f"✅ Mutation simulation log saved to {csv_filepath}")

# -------------------- AlphaFold JSON Generation --------------------
def generate_alphafold_batch_json(snapshot_dir, fasta_prefix, num_snapshots, output_json_file):
    job_list = []
    print(f"\n📦 Generating AlphaFold batch JSON from snapshots...")
    for i in range(num_snapshots):
        filename = os.path.join(snapshot_dir, f"{fasta_prefix}_{i:02d}.fasta")
        if not os.path.exists(filename):
            print(f"  ⚠️ Warning: File not found: {filename}. Skipping.")
            continue
        try:
            record = SeqIO.read(filename, "fasta")
            sequence = str(record.seq).upper()
            job_name = f"{fasta_prefix}_{i:02d}_EvoSim"
            job_list.append({
                "name": job_name,
                "modelSeeds": [],
                "sequences": [{"proteinChain": {"sequence": sequence, "count": 1}}]
            })
        except Exception as e:
            print(f"  ⚠️ Error processing {filename}: {e}")
    json_path = os.path.join(snapshot_dir, output_json_file)
    with open(json_path, 'w') as f:
        json.dump(job_list, f, indent=2)
    print(f"\n✅ AlphaFold batch JSON generated: {json_path}")

# -------------------- Main --------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="MAFFT → Ancestral → Mutation → AlphaFold JSON pipeline")
    parser.add_argument("input_fasta", help="Input multi-FASTA file with amino acid sequences")
    parser.add_argument("-m", "--mafft_path", default="mafft", help="Path to MAFFT executable")
    parser.add_argument("-a", "--mafft_options", default="--auto", help="MAFFT options")
    parser.add_argument("-o", "--aligned", default="aligned.fasta", help="Output aligned FASTA file")
    parser.add_argument("-c", "--consensus", default="Consensus_Sequence.fasta", help="Consensus sequence output file")
    parser.add_argument("--ancestor_id", default="Pigeon_HBB_Ancestor", help="ID for ancestral sequence record")
    parser.add_argument("--tree_file", default="calculated_tree.xml", help="Output PhyloXML tree")
    parser.add_argument("--newick_file", default="calculated_tree.nwk", help="Output Newick tree")
    parser.add_argument("--snapshot_dir", default="snapshots", help="Directory for mutation snapshots")
    parser.add_argument("--num_snapshots", type=int, default=20, help="Number of mutation snapshots")
    parser.add_argument("--mutation_rate", type=float, default=0.01, help="Mutation rate per site per snapshot")
    parser.add_argument("--lethal_rejection", type=float, default=0.95, help="Probability lethal changes are rejected")
    parser.add_argument("--json_file", default="alphafold_batch_submission.json", help="Output JSON for AlphaFold batch")

    args = parser.parse_args()

    aligned_file = run_mafft(args.input_fasta, args.aligned, args.mafft_path, args.mafft_options)
    reconstruct_weighted_ancestor(aligned_file, args.consensus, args.ancestor_id, args.tree_file, args.newick_file)
    simulate_mutations(args.consensus, args.snapshot_dir, args.num_snapshots, args.mutation_rate, args.lethal_rejection)
    generate_alphafold_batch_json(args.snapshot_dir, "protein_snap", args.num_snapshots, args.json_file)

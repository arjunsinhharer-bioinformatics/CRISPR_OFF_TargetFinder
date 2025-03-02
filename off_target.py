
"""
Created on Sun Mar 17 17:33:52 2024

@author: Arjunsinh Harer
"""

from Bio import SeqIO
from Bio.Seq import Seq
from pathlib import Path


def calculate_mismatches(seq1, seq2):
    """Calculate the number of mismatches between two sequences."""
    mismatches = sum(1 for a, b in zip(seq1, seq2) if a != b)
    return mismatches

def calculate_on_target_score(mismatches, length):
    """Calculate an on-target score. The fewer the mismatches, the higher the score.
       Simple scoring: 100 points, subtract 10 points for each mismatch."""
    score = max(100 - mismatches * 10, 0)
    return score

def read_gRNAs_from_file(file_path):
    """Read guide RNA sequences from a file, each gRNA on a new line."""
    with open(file_path, 'r') as file:
        gRNAs = [line.strip() for line in file if line.strip()]
    return gRNAs

def find_on_and_off_targets_for_gRNA_list(gRNAs, fasta_file, max_mismatches=3):
    """
    Find on-target sites within a FASTA file for a list of gRNAs, calculate on-target scores, 
    and identify potential off-targets with mismatches within the threshold. Checks both the
    gRNA sequence and its reverse complement.
    
    Parameters:
    - gRNAs: List of guide RNA sequences.
    - fasta_file: Path to the FASTA file containing target sequences.
    - max_mismatches: Maximum allowed mismatches for considering an off-target (int).
    
    Returns:
    - Dictionary with gRNA as key and a nested dictionary with 'on-target' and 'off-targets'.
    """
    targets = {gRNA: {'on-target': None, 'off-targets': []} for gRNA in gRNAs}
    
    for gRNA in gRNAs:
        gRNA_seq = Seq(gRNA)
        gRNA_seq_rc = gRNA_seq.reverse_complement()  # Reverse complement of the gRNA
        for record in SeqIO.parse(fasta_file, "fasta"):
            for i in range(len(record.seq) - len(gRNA_seq) + 1):
                segment = record.seq[i:i+len(gRNA_seq)]
                # Check both gRNA and its reverse complement against the segment
                for gRNA_variant in (gRNA_seq, gRNA_seq_rc):
                    mismatches = calculate_mismatches(gRNA_variant, segment)
                    score = calculate_on_target_score(mismatches, len(gRNA_variant))
                    if mismatches == 0:
                        targets[gRNA]['on-target'] = (str(segment), i, score)
                    elif mismatches <= max_mismatches:
                        targets[gRNA]['off-targets'].append((str(segment), i, mismatches, score))
    return targets



# Directory where the script is located
script_dir = Path(__file__).parent

# Find the first .txt file in the script directory that starts with "guides_"
guide_files = list(script_dir.glob("guides_*.txt"))
if guide_files:
    guides_file_path = guide_files[0]  # Take the first matching file
    gene_name = guides_file_path.stem.split("_")[1]  # Extract gene name

    # Construct the FASTA file path using the extracted gene name
    fasta_file = script_dir / f"{gene_name}_homo_sapiens_1.fasta"
else:
    raise FileNotFoundError("No guide RNA file found with the expected nomenclature ('guides_geneName.txt').")

# Now, guides_file_path and fasta_file can be used in your script as before


# Usage of guides_file_path and fasta_file in the rest of your script


# Read gRNAs from file
gRNAs = read_gRNAs_from_file(guides_file_path)

# Find targets for gRNAs read from the file
targets_dict = find_on_and_off_targets_for_gRNA_list(gRNAs, fasta_file)

for gRNA, target_data in targets_dict.items():
    print(f"Results for gRNA: {gRNA}")
    on_target = target_data['on-target']
    off_targets = target_data['off-targets']
    
    if on_target:
        sequence, position, score = on_target
        print(f"On-target sequence: {sequence}, Position: {position}, On-target score: {score}")
    else:
        print("No on-target match found.")
        
    if off_targets:
        print("Potential off-targets:")
        for off_target in off_targets:
            sequence, position, mismatches, score = off_target
            print(f"Off-target sequence: {sequence}, Position: {position}, Mismatches: {mismatches}, Score: {score}")
    else:
        print("No off-targets found.")
    print("-" * 50)



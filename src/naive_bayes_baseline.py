import pandas as pd
from Bio import SeqIO
from collections import Counter
import os
import time

# config
DATA_DIR = "./data/Train" 
OUTPUT_DIR = "./submissions"
NUM_PREDICTIONS = 1500


# LOAD AND ANALYSE TRAINING LABELS
def get_terms_frequencies(file_path):
    """Load the ground truth and counts the frequencies of each GO term."""
    print("1. Loading training terms...")

    # The file is TSV (tab separated) and has columns: EntryID, term, aspect
    # only care about the 'term' column for naive baseline
    train_terms = pd.read_csv(
        file_path,
        sep='\t',
        usecols=['EntryID', 'term'],
        header=0 # first row is header
    )

    # use counter for fast freq calc
    all_terms = train_terms['term'].tolist()
    term_counts = Counter(all_terms)

    # Get the most common terms and their freq/prob
    # Use the freq as the 'Score' (confidence)
    most_common = term_counts.most_common(NUM_PREDICTIONS)

    # Convert to a DataFrame for easier handling
    term_freq_df = pd.DataFrame(most_common, columns=['go_term', 'frequency'])
    total_annotations = term_freq_df['frequency'].sum()
    term_freq_df['score'] = term_freq_df['frequency'] / total_annotations

    print(f" -> Total unique GO terms: {len(term_counts):,}")
    print(f" -> using the top {NUM_PREDICTIONS} terms for predictions.")
    return term_freq_df.drop(columns=['frequency'])


# load test seq (testsuperset.fasta)
def get_test_proteins(file_path):
    """Parses the FASTA file to get a list of protein IDs."""
    print("2. Loading test proteins...")
    test_proteins = []
    
    # Biopython's SeqIO for FASTA parsing
    for record in SeqIO.parse(file_path, "fasta"):
        # "The ID in FASTA is formatted like "sp|ID|NAME"
        # We need he unique ID part
        parts = record.id.split('|')
        protein_id = parts[1] if len(parts) > 1 else record.id
        test_proteins.append(protein_id)
    
    print(f" -> total test proteins: {len(test_proteins):,}")
    return test_proteins

def generative_naive_submission(test_proteins, term_scores):
    """
    Predicts the same set of high-freq terms for every test protein.
    """
    print("3. Generating naive submission...")

    # Create the submission structure by combining every protein ID
    # with every predicted GO term/score.

    # Create a list of all protein rows
    submission_rows = []
    for protein_id in test_proteins:
        for inde, row in term_scores.iterrows():
            # protein_id, go_term, score (probability)
            submission_rows.append([
                protein_id,
                row['go_term'],
                f"{row['score']:.6f}" # format score to 6 dec places
            ])

    # Convert to DataFrame
    submission_df = pd.DataFrame(submission_rows)

    # Save the submission_file
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    timestamp = time.strftime("%Y%m%d-%H%M%S")
    output_file = os.path.join(OUTPUT_DIR, f"naive_baseline_{timestamp}.tsv")

    # submission file must be TSV, no index and header
    submission_df.to_csv(output_file, sep='\t', index=False, header=False)

    print(f"4. Success! Submission file saved to: {output_file}")
    print(f" -> Total prediction rows: {len(submission_df):, }")


if __name__ == "__main__":

    # File paths
    TRAIN_TERM_PATH = os.path.join(DATA_DIR, "train_terms.tsv")
    TEST_FASTA_PATH = os.path.join("./data/Test", "testsuperset.fasta")

    # 1. Get the naive scores
    naive_scores = get_terms_frequencies(TRAIN_TERM_PATH)

    # 2. Get the test protein IDs
    test_ids = get_test_proteins(TEST_FASTA_PATH)

    # 3. Create the sub file
    generative_naive_submission(test_ids, naive_scores)

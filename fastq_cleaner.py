import os
import gzip
from Bio import SeqIO
import matplotlib.pyplot as plt

# --- Part 1: Filtering and Counting ---
def analyze_fastq(filename):
    # dictionary to hold our unique sequences and their counts
    seq_counts = {}
    lengths = []
    n_percents = []

    # handle both zipped and unzipped files so we don't have to extract manually
    if filename.endswith(".gz"):
        handle = gzip.open(filename, "rt")
    else:
        handle = open(filename, "r")

    # read through the file using biopython
    for record in SeqIO.parse(handle, "fastq"):
        seq = str(record.seq).upper()
        seq_len = len(seq)
        
        # skip empty reads so the math doesn't crash later
        if seq_len == 0:
            continue
            
        # calculate the N percentage
        n_count = seq.count('N')
        n_percent = (n_count / seq_len) * 100
        
        # tool 2: throw out anything with more than 5% Ns
        # (check with saathvik if we need to change this threshold later)
        if n_percent > 5.0:
            continue 
            
        # tool 1: redundancy reduction
        # if the sequence is good, add it to our dictionary count
        seq_counts[seq] = seq_counts.get(seq, 0) + 1
            
        # save the stats for the graphs
        lengths.append(seq_len)
        n_percents.append(n_percent)

    # always close the file when done
    handle.close() 
    return seq_counts, lengths, n_percents

# --- Part 2: Sorting ---
def sort_by_abundance(seq_dict):
    # sorts the dictionary so the highest counts are at the top
    return sorted(seq_dict.items(), key=lambda x: x[1], reverse=True)

# --- Part 3: Saving Files and Plotting ---
def save_reports(sorted_seqs, text_file, fasta_file, total_reads):
    # save the text dashboard first
    with open(text_file, 'w') as t_file:
        t_file.write("--- QUALITY DASHBOARD & ABUNDANCE TABLE ---\n")
        t_file.write(f"Total reads passing 5% filter: {total_reads}\n")
        t_file.write(f"Unique sequences found: {len(sorted_seqs)}\n\n")
        
        # only writing the top 50 so the file isn't massive
        for i, (seq, count) in enumerate(sorted_seqs[:50]): 
            t_file.write(f"Rank {i+1} | Count: {count} | Seq: {seq}\n")

    # save the fasta file for downstream analysis
    with open(fasta_file, 'w') as f_file:
        for i, (seq, count) in enumerate(sorted_seqs):
            # standard fasta format needs the > at the start
            f_file.write(f">unique_seq_{i+1}_abundance_{count}\n")
            f_file.write(f"{seq}\n")

def create_plots(lengths, n_percents, sorted_seqs, folder_name):
    # set up a wide window for all 3 charts
    plt.figure(figsize=(15, 6)) 
    
    # main title for the whole thing
    plt.suptitle(f"Analysis Results for: {folder_name}", fontsize=16, fontweight='bold')
    
    # Chart 1: Bar chart for top 10 sequences
    plt.subplot(1, 3, 1)
    top_10 = sorted_seqs[:10]
    
    # make simple labels like Seq 1, Seq 2 so it fits on the axis
    x_labels = [f"Seq {i+1}" for i in range(len(top_10))]
    counts = [item[1] for item in top_10]
    
    plt.bar(x_labels, counts, color='lightgreen', edgecolor='black')
    plt.title(f'Top 10 Abundance ({folder_name})')
    plt.xticks(rotation=45) # tilt labels
    plt.ylabel('Read Count')
    
    # Chart 2: Read length histogram
    plt.subplot(1, 3, 2)
    plt.hist(lengths, bins=30, color='skyblue', edgecolor='black')
    plt.title(f'Read Lengths ({folder_name})')
    plt.xlabel('Length (bp)')
    plt.ylabel('Frequency')
    
    # Chart 3: N-content histogram
    plt.subplot(1, 3, 3)
    plt.hist(n_percents, bins=30, color='salmon', edgecolor='black')
    plt.title(f'N-Content <5% ({folder_name})')
    plt.xlabel('% of N bases')
    
    # fix the spacing so the titles don't overlap
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()

# --- Main Program ---
print("--- Preprocessing Tools ---")

# ask which folder to run
target_folder = input("Enter the folder name to analyze (e.g., barcode09): ").strip()

# make sure the folder actually exists before trying to run
if os.path.exists(target_folder) and os.path.isdir(target_folder):
    
    # naming our output files based on the folder we typed
    dashboard_file = f"{target_folder}_dashboard.txt"
    fasta_out = f"{target_folder}_deduplicated.fasta"
    
    # variables to hold data across all files in the folder
    final_counts = {}
    final_lengths = []
    final_n_percents = []
    total_good_reads = 0

    print(f"\nProcessing files inside '{target_folder}'...")

    # loop through everything in the directory
    for file_name in os.listdir(target_folder):
        # only grab fastq or compressed fastq files
        if file_name.endswith((".fastq", ".fq", ".gz")):
            file_path = os.path.join(target_folder, file_name)
            print(f" > Reading: {file_name}")
            
            # run our function
            counts, lens, ns = analyze_fastq(file_path)
            
            # add the results to our totals
            total_good_reads += len(lens)
            final_lengths.extend(lens)
            final_n_percents.extend(ns)
            
            # merge the dictionaries
            for s, c in counts.items():
                final_counts[s] = final_counts.get(s, 0) + c
    
    # sort the final combined data
    sorted_data = sort_by_abundance(final_counts)
    
    # save everything and show the graphs
    save_reports(sorted_data, dashboard_file, fasta_out, total_good_reads)
    
    print(f"\nDone! Total good reads: {total_good_reads}")
    print(f"Saved dashboard to: {dashboard_file}")
    print(f"Saved FASTA to: {fasta_out}")
    print("Opening the charts...")
    
    create_plots(final_lengths, final_n_percents, sorted_data, target_folder)

else:
    print(f"\nError: Couldn't find a folder named '{target_folder}'. Check spelling!")
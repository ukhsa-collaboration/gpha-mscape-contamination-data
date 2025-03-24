#!/usr/bin/env python

import os
import pandas as pd
import numpy as np

spikeins = {
        12242:["Tobamovirus","Tobacco_mosaic_virus"]
    }

def convert_to_numeric(column):
        if column.name not in ['Rank', 'Scientific_Name', "Domain"]:
            return pd.to_numeric(column, errors='coerce')
        else:
            return column

def save_perc_and_count_dfs(needed_samples, reports):
    # Initialize an empty list to store dataframes
    dfs_perc = []
    dfs_count = []

    #   Loop over each kraken file in reports directory
    for sample in needed_samples:
        found = False
        for filename in reports:
            if f'{sample}' in filename:
                found = True

                #open new file and read it line by line
                file = open(filename)

                #create 3 lists for count of sequences then rank and scientific name
                perc_seqs = []
                read_counts = []
                rank = []
                sci_name = []
                domain = []

                #only start reading each file when it starts listing our domain
                current_domain = None

                # Iterate over each line in the file
                for line in file:
                    if line.startswith("%"):
                        continue

                    # Split the line into columns
                    columns = line.split()

                    #record the current scientific name of this line
                    current_name = columns[5]
                    current_rank = columns[3]

                    # Check for the microbial category in scientific name column and set the flag to start readin
                    if current_name in ["Sar", "Discoba", "Metamonada"]:
                        current_domain = "Protists"
                    elif current_name == "Fungi":
                        current_domain = "Fungi"
                    elif current_name == "Bacteria":
                        current_domain = "Bacteria"
                    elif current_name == "Viruses":
                        current_domain = "Viruses"
                    elif current_name == "Archaea":
                        current_domain = "Archaea"
                    elif (current_domain == "Protists" and current_rank == "D1") \
                        or (current_domain == "Fungi" and current_rank == "K") \
                        or (current_domain in ["Bacteria","Viruses","Archaea"] and current_rank == "D"):
                        current_domain = None

                    perc_seqs.append(columns[0])
                    read_counts.append(columns[1])
                    rank.append(current_rank)
                    sci_name.append('_'.join(columns[5:])) #this turns scientific name (which sometimes have multiple words) into a list within a list
                    domain.append(current_domain)

                # Extract the sample ID from sample (climb_id)
                sample_ID = sample

                # Turn the four lists into a dataframe, using sample ID in place of "% of seqs" or "read counts", depending on whether you want counts or percentages
                df_perc = pd.DataFrame({sample_ID: read_counts, "Rank": rank, "Scientific_Name": sci_name, "Domain":domain})
                df_count = pd.DataFrame({sample_ID: read_counts, "Rank": rank, "Scientific_Name": sci_name, "Domain":domain})

                # Set the index to scientific name for merging later
                df_perc.set_index("Scientific_Name", inplace=True)
                df_count.set_index("Scientific_Name", inplace=True)

                #add the new dataframe to the list of dataframes
                dfs_perc.append(df_perc)
                dfs_count.append(df_count)


    # Merge the DataFrames on a specific column
    merged_perc_df = pd.concat(dfs_perc, axis = 0, join= "outer")  # Change join to 'outer' for outer join
    merged_perc_df = merged_perc_df.apply(convert_to_numeric).fillna(0)

    merged_count_df = pd.concat(dfs_count, axis = 0, join= "outer")  # Change join to 'outer' for outer join
    merged_count_df = merged_count_df.apply(convert_to_numeric).fillna(0)

    # Remove spikeins
    for spike in spikeins:
        merged_perc_df = merged_perc_df.loc[~merged_perc_df.index.isin(spikeins[spike])]
        merged_count_df = merged_count_df.loc[~merged_count_df.index.isin(spikeins[spike])]

    # Save to file
    merged_perc_df.to_csv('perc.csv', index=True)
    merged_count_df.to_csv('counts.csv', index=True)

def main():
    needed_samples = ["NTC_231219","Ntc_231122","Ntc_240105"]
    reports = ["/Users/rmcolq/Downloads/neg_control_reports/reports/GSTT_NTC/kraken2/NTC_231219.kraken2.report.txt","/Users/rmcolq/Downloads/neg_control_reports/reports/GSTT_NTC/kraken2/Ntc_231213.kraken2.report.txt","/Users/rmcolq/Downloads/neg_control_reports/reports/GSTT_NTC/kraken2/NTC_240208.kraken2.report.txt","/Users/rmcolq/Downloads/neg_control_reports/reports/GSTT_NTC/kraken2/Ntc_240105.kraken2.report.txt","/Users/rmcolq/Downloads/neg_control_reports/reports/GSTT_NTC/kraken2/Ntc_231122.kraken2.report.txt"]
    save_perc_and_count_dfs(needed_samples, reports)

if __name__ == "__main__":
    main()
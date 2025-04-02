#!/usr/bin/env python
import pandas as pd
import os

spikeins = {
        12242:["Tobamovirus","Tobacco_mosaic_virus"]
    }

def get_label(ids, site_key, run_id=None):
    ids_list = [id.lower() for id in ids]

    ids_list_anon = []
    for id in ids_list:
        if id in site_key:
            ids_list_anon.append(site_key[id])
        else:
            ids_list_anon.append(id)
    ids_list = ids_list_anon

    ids_list = ' '.join(ids_list)
    if "public" in ids and run_id:
        ids_list = ids_list.replace('public', 'public_'+run_id) #add run_id to name
    ids_list = ids_list.replace('water_extraction_control', '(water)')
    ids_list = ids_list.replace('resp_matrix_mc110', '(matrix)')
    return ids_list

def convert_to_numeric(column):
        if column.name not in ['Rank', 'Scientific_Name']:
            return pd.to_numeric(column, errors='coerce')
        else:
            return column

#Create a dataframe for average counts per dataset in each microbe type(taxon) in "bacteria", "fungi", "viruses", "archaea", and "protists"
def save_perc_and_count_dfs(needed_samples, reports, df_type):
    # Initialize an empty list to store dataframes
    merged_df = None

    #Loop over each kraken file in reports directory
    for sample in needed_samples:
        found = False
        for filename in os.listdir(reports):
            if f'{sample}' in filename:
                found = True
                #open new file and read it line by line
                file = open(f'{reports}/{filename}')

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

                # Turn the four lists into a dataframe, using sample ID in place of "% of seqs" or "read counts", depending on whether you want counts or percentages
                if df_type == "perc":
                    df = pd.DataFrame({sample: perc_seqs, "Rank": rank, "Scientific_Name": sci_name, "Domain":domain})
                else:
                    df = pd.DataFrame({sample: read_counts, "Rank": rank, "Scientific_Name": sci_name, "Domain":domain})

                # Set the index to scientific name for merging later
                df.set_index(["Scientific_Name", "Rank", "Domain"], inplace=True)
                #add the new dataframe to the merged df
                if merged_df is None:
                    merged_df = df
                else:
                    merged_df = merged_df.join(df, how="outer")
        if not found:
            print(f"No kraken report for sample {sample} has been provided!")


    # Merge the DataFrames on a specific column
    merged_df = merged_df.apply(convert_to_numeric).fillna(0)
    merged_df = merged_df.reset_index()

    return merged_df

    # Remove spikeins
    #for spike in spikeins:
    #    merged_perc_df = merged_perc_df.loc[~merged_perc_df.index.isin(spikeins[spike])]
    #    merged_count_df = merged_count_df.loc[~merged_count_df.index.isin(spikeins[spike])]

    # Save to file
    #merged_perc_df.to_csv('perc.csv', index=True)
    #merged_count_df.to_csv('counts.csv', index=True)
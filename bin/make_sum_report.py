#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
import io
import base64
from mako.template import Template
import numpy as np
import seaborn as sns
import os
import argparse
import json

from utils import spikeins, convert_to_numeric, get_label

#Create a dataframe for average counts per dataset in each microbe type(taxon) in "bacteria", "fungi", "viruses", "archaea", and "protists"
def get_each_taxon(needed_samples, reports, microbe_type):
    # Initialize an empty list to store dataframes
    dfs = []

    #Loop over each kraken file in reports directory
    for sample in needed_samples:
        found = False
        for filename in reports:
            if f'{sample}' in filename:
                found = True
                #open new file and read it line by line
                file = open(filename)
            
                #create 3 lists for count of sequences then rank and scientific name
                read_counts = []
                rank = []
                sci_name = []   
        
                # Iterate over each line in the file
                for line in file:
                    if line.startswith("%"):
                        continue

                    # Split the line into columns
                    columns = line.split()

                    #record the current scientific name of this line
                    current_name = columns[5]
                    current_id = int(columns[4])
                        
                    # Check for the microbial category in scientific name column and set the flag to start readin
                    if microbe_type == "Protists":
                        if current_name in ["Sar", "Discoba", "Metamonada"]:
                            read_counts.append(line.split()[1])
                            rank.append(line.split()[3])            
                        #this turns scientific name (which sometimes have multiple words) into a list within a list
                            sci_name.append(line.split()[5:])
                    elif microbe_type == "All":
                        if current_name == "root" or current_id in spikeins:
                            read_counts.append(line.split()[1])
                            rank.append(line.split()[3])            
                        #this turns scientific name (which sometimes have multiple words) into a list within a list
                            sci_name.append(line.split()[5:])
                    # Get list of taxa in bacteria, virus, archaea domain 
                    elif microbe_type == "Viruses":
                        if current_name == microbe_type or current_id in spikeins:
                            read_counts.append(line.split()[1])
                            rank.append(line.split()[3])            
                        #this turns scientific name (which sometimes have multiple words) into a list within a list
                            sci_name.append(line.split()[5:])
                            #print("Appended line in kraken report:")
                            #print(line)
                    else:
                        if current_name == microbe_type:
                            read_counts.append(line.split()[1])
                            rank.append(line.split()[3])            
                        #this turns scientific name (which sometimes have multiple words) into a list within a list
                            sci_name.append(line.split()[5:])
                                            
                # Extract the sample ID from sample (climb_id)
                sample_ID = sample
            
                # Turn the four lists into a dataframe, using sample ID in place of "% of seqs" or "read counts", depending on whether you want counts or percentages
                df_new_file = pd.DataFrame({sample_ID: read_counts, "Rank": rank, "Scientific_Name": sci_name})
                
                #add the new dataframe to the list of dataframes
                dfs.append(df_new_file)

        if not found:
            print(f"No kraken report for sample {sample} has been provided!")
    
    
    # Merge the DataFrames on a specific column
    merged_df = pd.concat(dfs, axis = 0, join= "outer")  # Change join to 'outer' for outer join
    #turn scientific_name from a list to a string
    merged_df["Scientific_Name"] = ['_'.join(lst) for lst in merged_df["Scientific_Name"]]
    #merge on scientific name - this means that no scientific name is repeated if it's present in different dataframes
    merged_on_sci = merged_df.groupby('Scientific_Name', as_index=False).first()

    # Apply the function to each column
    numeric_df = merged_on_sci.apply(convert_to_numeric)
    
    #change NaN to "0"
    numeric_df = numeric_df.fillna(0)

    #Get list of spikeins and make a spikeins dataframe from main dataframe
    spikein_names = []
    for spike in spikeins:
        spike_list = spikeins[spike]
        spikein_names.extend(spike_list)
    
    if isinstance(spikein_names, str):
        spike_df = numeric_df.loc[numeric_df["Scientific_Name"] == spikein_names]
    else:
        spike_df = numeric_df.loc[numeric_df["Scientific_Name"].isin(spikein_names)]
    
    #if microbe_type == "Viruses":
    #    print("Dataframe with all potential spikeins:")
    #    print(spike_df)

    #Get sum of spikeins
    spike_df = spike_df.drop(columns=['Scientific_Name', 'Rank'])
    copy_df = spike_df.copy()
    copy_df["Sum"] = copy_df.sum(axis=1) #sum for each spikein taxon
    spike_sum = copy_df["Sum"].sum(axis=0) #sum for all spikein taxa combined
     
    #Filtering the dataframe to find the total of whichever taxon    
    if microbe_type != "Protists" and microbe_type != "All":
        taxon_needed = numeric_df.loc[numeric_df["Scientific_Name"] == microbe_type]
    elif microbe_type == "All":
        taxon_needed = numeric_df.loc[numeric_df["Scientific_Name"] == "root"]
    else:
        taxon_needed = numeric_df.loc[numeric_df["Scientific_Name"].isin(["Sar", "Discoba", "Metamonada"])]

    #if microbe_type == "Viruses":
    #    print("Dataframe for domain level Virus counts:")
    #    print(taxon_needed)

    # Select columns to sum (excluding 'Kingdom' and 'Scientific Name')
    columns_needed = taxon_needed.drop(columns=['Scientific_Name', 'Rank'])
    columns_copy = columns_needed.copy()
    #Get the total counts in each dataset for the phylum we are checking for 
    columns_sum = columns_copy.sum(axis=1)
    #Get the sum of all taxa that we are checking for - needed for protists as there are 3 taxa
    columns_copy['sum'] = columns_sum
    microbe_sum = columns_copy['sum'].sum(axis=0)

    real_count = int(microbe_sum) - int(spike_sum)

    if microbe_type == "Viruses":
        print(f'Count with spikeins removed: all virus counts({microbe_sum}) - all tobamo virus counts ({spike_sum}) = {real_count}')
    
    #Get the total number of columns/samples
    columns_no = columns_copy.shape[1]

    #Get average count per sample in a dataset
    average = real_count/columns_no

    data = {'Name': microbe_type, 'Average Count': average}
    current_df = pd.DataFrame(data, index=[0])
    return current_df

#Merge all taxon count dataframes together to get all microbe type count data per dataset
def get_microbial_load(set, needed_samples, reports):
    taxa = ["All", "Bacteria", "Fungi", "Viruses", "Archaea", "Protists"]
    current_dfs = []
    for microbe_type in taxa:
        current_df = get_each_taxon(needed_samples, reports, microbe_type)
        current_dfs.append(current_df)
    # Merge the DataFrames on a specific column
    merged_df = pd.concat(current_dfs, axis = 0, join= "outer")  # Change join to 'outer' for outer join
    #merge on scientific name - this means that no scientific name is repeated 
    merged_df = merged_df.groupby('Name', as_index=True).first()
    transposed_df = merged_df.transpose()
    transposed_df.reset_index(inplace=True)

    transposed_df['index'] = set
    return transposed_df

#Make microbe count table for each dataset and merge them together
def make_microbial_count_table(reports, grouped_metadata, site_key):

    datasets = [] #total datasets
    samples = [] #samples grouped by datasets

    for sets in grouped_metadata:
            ids = list(sets[0])
            if "public" in ids[0].lower(): #if dataset is public
                subgroup = sets[1].groupby(['run_id'])
                for subset in subgroup:
                    run_id = ''.join(subset[0])
                    #turn scientific_name from a list to a string
                    ids_list = get_label(ids, site_key, run_id=run_id)

                    datasets.append(ids_list)
                    table = subset[1] #list of all ids in sub_dataset
                    samples.append(list(table['climb_id'])) #climb id
            else:
                #turn scientific_name from a list to a string
                ids_list = get_label(ids, site_key)
                datasets.append(ids_list)
                table = sets[1] #list of all ids in dataset
                samples.append(list(table['climb_id'])) #climb id

    loop = 0
    single_dfs = []
    for set in datasets:
        needed_samples = samples[loop] #our current set of sample names
        transposed_df = get_microbial_load(set, needed_samples, reports) #get dataset name, sample names, and all reports
        single_dfs.append(transposed_df)
        loop += 1
        
    #merge all dataframes together
    final_table = pd.concat(single_dfs, axis = 0, join = "outer")
    final_table = final_table.groupby('index', as_index=False).first()
    final_table.fillna(value=0.0, inplace=True)
    return final_table

def get_species_count(needed_samples, reports, microbe_type, taxon_level, filter_count):
    # Initialize an empty list to store dataframes
    dfs = []
#   Loop over each kraken file in reports directory
    for sample in needed_samples:
        for filename in reports:
            if f'{sample}' in filename:
                #open new file and read it line by line
                file = open(filename)
        
                #create 3 lists for count of sequences then rank and scientific name
                read_counts = []
                rank = []
                sci_name = []
            
            #only start reading each file when it starts listing our domain
                start_reading = False    

                # Iterate over each line in the file
                for line in file:
                    # Split the line into columns
                    columns = line.split()
        
                    #record the current scientific name of this line
                    current_name = columns[5]
                    current_rank = columns[3]
    
                    # Check for the microbial category in scientific name column and set the flag to start readin
                    if microbe_type == "Protists":
                        if current_name in ["Sar", "Discoba", "Metamonada"]:
                            start_reading = True
                        else:
                            if current_rank == "D1":
                                start_reading = False
                    # Get list of taxa in Fungi kingdom only
                    elif microbe_type == "Fungi":
                        if current_name == microbe_type:
                            start_reading = True
                        else:
                            if current_rank == "K":
                                start_reading = False
                    # Get list of taxa in bacteria, virus, archaea domain
                    elif microbe_type == f'Bacteria > {filter_count}':
                        if current_name == "Bacteria":
                            start_reading = True
                        else:
                            if current_rank == "D":
                                start_reading = False  
                    else:
                        if current_name == microbe_type:
                            start_reading = True
                        else:
                            if current_rank == "D":
                                start_reading = False
                        
                                    
                    if start_reading:
                        read_counts.append(line.split()[1])
                        rank.append(line.split()[3])            
                    #this turns scientific name (which sometimes have multiple words) into a list within a list
                        sci_name.append(line.split()[5:])

                # Extract the sample ID from sample (climb_id)                           
                sample_ID = sample
            
                # Turn the four lists into a dataframe, using sample ID in place of "% of seqs" or "read counts", depending on whether you want counts or percentages
                df_new_file = pd.DataFrame({sample_ID: read_counts, "Rank": rank, "Scientific_Name": sci_name})
                
                #add the new dataframe to the list of dataframes
                dfs.append(df_new_file)


    # Merge the DataFrames on a specific column
    merged_df = pd.concat(dfs, axis = 0, join= "outer")  # Change join to 'outer' for outer join
    #turn scientific_name from a list to a string
    merged_df["Scientific_Name"] = ['_'.join(lst) for lst in merged_df["Scientific_Name"]]
    #merge on scientific name - this means that no scientific name is repeated if it's present in both the gut and lung dataframes
    merged_on_sci = merged_df.groupby('Scientific_Name', as_index=False).first()

    # Apply the function to each column
    numeric_df = merged_on_sci.apply(convert_to_numeric)   
    #change NaN to "0"
    numeric_df = numeric_df.fillna(0)
    
    #remove spikeins
    for spike in spikeins:
        numeric_df = numeric_df.loc[~numeric_df["Scientific_Name"].astype(str).isin(spikeins[spike])]

    # Define keywords and columns to search
    keywords = [taxon_level] #the taxonomy level I want to filter for
    columns_to_search = ['Rank']
    
    # Boolean indexing to filter rows showing only genus in the rank column
    filtered_df = numeric_df[numeric_df[columns_to_search].apply(lambda x: x.isin(keywords).any(), axis=1)]

    if len(filtered_df.index) > 0 :
        # Rearrange column 'Scientific Name' to the first position
        wordy_columns = ['Scientific_Name', 'Rank']
        # check if columns are not in the wordy_columns list
        column_order = ['Scientific_Name'] + ['Rank'] + [col for col in filtered_df.columns if col not in wordy_columns]
        filtered_df = filtered_df[column_order]

        if microbe_type == f"Bacteria > {filter_count}":
            # Select columns to sum (excluding 'Kingdom' and 'Scientific Name')
            columns_to_sum = filtered_df.drop(columns=['Scientific_Name', 'Rank'])            
            # Calculate the sum of values in each row
            column_sums = columns_to_sum.sum(axis=1)
            # Count the number of columns
            num_columns = columns_to_sum.shape[1]
            #either get a sum of all seqs, or normalise it by dividing by num_columns
            filtered_df["Average_Counts"] = column_sums/num_columns
            above_count = filtered_df[filtered_df["Average_Counts"] >= filter_count]

            number = above_count.shape[0]
        else:
            number = filtered_df.shape[0]
            
        data = {'Name': microbe_type, 'Number of Genera': number}
        df = pd.DataFrame(data, index=[0])
        return df

def get_all_taxa(set, needed_samples, reports, taxon_level, filter_count):
    taxa = ["Bacteria", f"Bacteria > {filter_count}", "Fungi", "Viruses", "Archaea", "Protists"]
    current_dfs = []
    for microbe_type in taxa:
        df = get_species_count(needed_samples, reports, microbe_type, taxon_level, filter_count)
        current_dfs.append(df)
        
    #merge all dataframes together
    table = pd.concat(current_dfs, axis = 0, join = "outer")
    table = table.groupby('Name', as_index=True).first()
    table.fillna(value=0.0, inplace=True)

    transposed_df = table.transpose()
    transposed_df.reset_index(inplace=True)

    transposed_df['index'] = set
    return transposed_df

def make_richness_table(reports, grouped_metadata, taxon_level, filter_count, site_key):
 
    datasets = [] #total datasets
    samples = [] #samples grouped by datasets

    for sets in grouped_metadata:
            ids = list(sets[0])
            if "public" in ids[0].lower(): #if dataset is public
                subgroup = sets[1].groupby(['run_id'])
                for subset in subgroup:
                    run_id = ''.join(subset[0])
                    #turn scientific_name from a list to a string
                    ids_list = get_label(ids, site_key, run_id=run_id)

                    datasets.append(ids_list)
                    table = subset[1] #list of all ids in sub_dataset
                    samples.append(list(table['climb_id'])) #climb id
            else:
                #turn scientific_name from a list to a string
                ids_list = get_label(ids, site_key)
                datasets.append(ids_list)
                table = sets[1] #list of all ids in dataset
                samples.append(list(table['climb_id'])) #climb id

    loop = 0
    single_dfs = []
    for set in datasets:
        needed_samples = samples[loop] #our current set of sample names
        transposed_df = get_all_taxa(set, needed_samples, reports, taxon_level, filter_count)
        single_dfs.append(transposed_df)
        loop += 1
        
    #merge all dataframes together
    final_table = pd.concat(single_dfs, axis = 0, join = "outer")
    final_table = final_table.groupby('index', as_index=False).first()
    final_table.fillna(value=0.0, inplace=True)

    return final_table 

def make_perc_df(needed_samples, reports):
    # Initialize an empty list to store dataframes
    dfs = []
#   Loop over each kraken file in reports directory
    for sample in needed_samples:
        for filename in reports:
            if f'{sample}' in filename:
                #open new file and read it line by line
                file = open(filename)
        
                #create 5 lists for %/count of sequences then rank and scientific name
                perc_seqs = []
                read_count = []
                rank = []
                sci_name = []
                
                # Iterate over each line in the file
                for line in file:
                    #make the lists
                    perc_seqs.append(line.split()[0])
                    read_count.append(line.split()[1])
                    rank.append(line.split()[3])                                
                    #this turns scientific name (which sometimes have multiple words) into a list within a list
                    sci_name.append(line.split()[5:])               
                                

                #Get total read count from "root"
                total_reads = read_count[1]

                # Turn the four lists into a dataframe, using sample ID in place of "% of seqs" or "read counts", depending on whether you want counts or percentages
                df_new_file = pd.DataFrame({total_reads: perc_seqs, "Rank": rank, "Scientific_Name": sci_name})
                
                #add the new dataframe to the list of dataframes
                dfs.append(df_new_file)
    
    # Merge the DataFrames on a specific column
    merged_df = pd.concat(dfs, axis = 0, join= "outer")  # Change join to 'outer' for outer join
    #turn scientific_name from a list to a string
    merged_df["Scientific_Name"] = ['_'.join(lst) for lst in merged_df["Scientific_Name"]]
    #merge on scientific name - this means that no scientific name is repeated if it's present in both the gut and lung dataframes
    merged_on_sci = merged_df.groupby('Scientific_Name', as_index=False).first()

    # Apply the function to each column
    numeric_df = merged_on_sci.apply(convert_to_numeric)
    
    # Select columns to sum (excluding 'Rank' and 'Scientific Name')
    columns_to_sum = numeric_df.drop(columns=['Rank', 'Scientific_Name'])
    
    # Calculate the sum of values in each row
    column_sums = columns_to_sum.sum(axis=1)
    # Count the number of columns
    num_columns = numeric_df.shape[1]
    #either get a sum of all seqs, or normalise it by dividing by num_columns
    numeric_df["Perc_Seqs_Overall"] = column_sums/num_columns
    
    #change NaN to "0"
    numeric_df = numeric_df.fillna(0)
    
    #remove spikeins
    for spike in spikeins:
        numeric_df = numeric_df.loc[~numeric_df["Scientific_Name"].astype(str).isin(spikeins[spike])]

    #Filtering the dataframe    
    # Define keywords and columns to search
    keywords = ['G'] #the taxonomy level I want to filter for
    columns_to_search = ['Rank']
    
    # Boolean indexing to filter rows showing only genus in the rank column
    genus_df = numeric_df[numeric_df[columns_to_search].apply(lambda x: x.isin(keywords).any(), axis=1)]
    
    # and remove rows with zero values (0 or 0.0)
    keywords_to_remove = [0]
    # Boolean indexing to filter rows in % of seqs overall that contain 0
    filtered_df = genus_df[~genus_df["Perc_Seqs_Overall"].astype(float).isin(keywords_to_remove)]
    
    # Rearrange column 'Scientific Name' to the first position
    wordy_columns = ['Scientific_Name']
    # check if columns are not in the wordy_columns list
    column_order = ['Scientific_Name'] + [col for col in filtered_df.columns if col not in wordy_columns]
    filtered_df = filtered_df[column_order]
    
    #drop column “rank” since they are all Fs, and turn each count into an integer
    no_rank_df = filtered_df.drop(columns=['Rank'])
    return no_rank_df

def get_heatmap(reports, grouped_metadata, site_key):
    
    #group by site
    datasets = []
    public_datasets = []
    samples = []
    public_samples = []

    for sets in grouped_metadata:
            ids = list(sets[0])
            if "public" in ids[0].lower(): #if dataset is public
                subgroup = sets[1].groupby(['run_id'])
                for subset in subgroup:
                    run_id = ''.join(subset[0])
                    #turn scientific_name from a list to a string
                    ids_list = get_label(ids, site_key, run_id=run_id)

                    datasets.append(ids_list)
                    public_datasets.append(ids_list)
                    table = subset[1] #list of all ids in sub_dataset
                    samples.append(list(table['climb_id'])) #climb id
                    public_samples.append(list(table['climb_id']))
            else:
                #turn scientific_name from a list to a string
                ids_list = get_label(ids, site_key)
                datasets.append(ids_list)
                table = sets[1] #list of all ids in dataset
                samples.append(list(table['climb_id'])) #climb id
            
    
    average_list = []
    
    loop = 0
    for set in datasets:
        needed_samples = samples[loop] #our current set of sample names
        perc_df = make_perc_df(needed_samples, reports)

        perc_df = perc_df.rename(columns={c: f"{set}_"+c for c in perc_df.columns if c not in ['Scientific_Name', 'Perc_Seqs_Overall']})

        #Change all "average percentage" columns to their respective dataset names to avoid clashes when merging
        perc_df[set] = perc_df["Perc_Seqs_Overall"]
        perc_df = perc_df.drop(columns=["Perc_Seqs_Overall"])
        average_list.append(perc_df)

        loop += 1

    #First define the merged dataframe by the starting dataframe
    loop_count = 0
    merge_df = average_list[0]
    #Then merge all other dataframes to the pre-existing dataframe
    for df in average_list:
        loop_count += 1
        if loop_count < len(datasets):
            current_df = average_list[loop_count]
            merge_df = merge_df.merge(current_df, on="Scientific_Name", how="outer")
    
    merge_df.fillna(value=0, inplace=True)

    mscape_datasets = []
    mscape_samples = []
    for set in datasets:
        if set not in public_datasets:
            mscape_datasets.append(set)
            mscape_samples.append(set)


    def sorting(merge_df, type_of_directories, dataset_names, all_dataset_names):
        # Select columns to sum ('Scientific Name')
        columns_to_sum = merge_df[dataset_names]

        # Calculate the sum of values in each row
        column_sums = columns_to_sum.sum(axis=1)
        merge_df['Average'] = column_sums/len(type_of_directories)

        # Calculate the max of values in each row
        # merge_df['Max'] = columns_to_sum.max(axis=1)

        sorted_df = merge_df.sort_values(by="Average", ascending=False)
        top_df = sorted_df.head(20)
    
        total_df = top_df.drop(columns=all_dataset_names)
        total_df = total_df.drop(columns=["Average"])

        both_dfs = []
        both_counts = []

        mscape_dfs = []
        publics = "public"

        #make specific dfs for all dataframes in mscape_datasets
        for dataset_name in mscape_datasets:
            name_col = "Scientific_Name"
            df = total_df.loc[:, total_df.columns.str.contains(dataset_name)]
            df[name_col] = total_df[name_col]
            mscape_dfs.append(df)
        
        loop_count = 0
        sorted_mscapes = []
        sorted_m_counts = []
        for mscape_df in mscape_dfs:
            #Get rid of dataset as prefixes
            mscape_df = mscape_df.rename(columns={c: c.replace(f"{mscape_datasets[loop_count]}_", "") for c in mscape_df.columns if c not in ['Scientific_Name']})
            mscape_df.sort_index(axis=1, inplace=True)
            sorted_mscapes.append(mscape_df)

            mscape_matrix = mscape_df.drop(columns=["Scientific_Name"])
            mscape_counts = mscape_matrix.transpose()
            mscape_counts = mscape_counts.reset_index()
            sorted_m_counts.append(mscape_counts['index'])
            loop_count += 1
        
        public_df = total_df.loc[:, total_df.columns.str.contains(publics, case=False)]
        public_df["Scientific_Name"] = total_df["Scientific_Name"]
        public_df = public_df.rename(columns={c: c.split("_")[-1] for c in public_df.columns if c not in ['Scientific_Name']})
        public_df.sort_index(axis=1, inplace=True)

        public_matrix = public_df.drop(columns=["Scientific_Name"])
        public_counts = public_matrix.transpose()
        public_counts = public_counts.reset_index()

        both_dfs.append(sorted_mscapes)
        both_dfs.append(public_df)

        both_counts.append(sorted_m_counts)
        both_counts.append(public_counts['index'])

        return both_dfs, both_counts

    sort_by_average, both_counts = sorting(merge_df, samples, datasets, datasets)
    sort_by_mscape, both_counts = sorting(merge_df, mscape_samples, mscape_datasets, datasets)

    return sort_by_average, sort_by_mscape, mscape_datasets, both_counts 

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse text files and a CSV file.")
    parser.add_argument('--reports', nargs='+', help="List of text files", required=True)
    parser.add_argument('--metadata', help="CSV file path", required=True)
    parser.add_argument('--site_key', help="JSON file specifying site name to number", required=True)
    parser.add_argument('--plots_dir', help="Shannon plots directory", required=True)
    parser.add_argument('--final_report', help="Output directory", required=True)
    parser.add_argument('--template', help="HTMl template", required=True)
    args = parser.parse_args()

    reports = args.reports
    metadata = pd.read_csv(args.metadata)
    plots_dir = args.plots_dir # Import R plots
    palette = pd.read_csv(f'{plots_dir}/colour_palette.txt', delimiter='\t')
    sns.set_style("whitegrid")
    all_directories = []
    mscape_directories = []
    #group by site and other identifiers
    grouped_metadata = metadata.groupby(['site', 'control_type_details'])
    #site name to number
    site_key = {}
    with open(args.site_key, 'r') as f:
        site_key = json.load(f)

    #group by site
    datasets = []
    public_datasets = []
    samples = []
    public_samples = []

    for sets in grouped_metadata:
            ids = list(sets[0])
            if "public" in ids[0].lower(): #if dataset is public
                subgroup = sets[1].groupby(['run_id'])
                for subset in subgroup:
                    run_id = ''.join(subset[0])
                    #turn scientific_name from a list to a string
                    ids_list = get_label(ids, site_key, run_id=run_id)

                    datasets.append(ids_list)
                    public_datasets.append(ids_list)
                    table = subset[1] #list of all ids in sub_dataset
                    samples.append(list(table['climb_id'])) #climb id
                    public_samples.append(list(table['climb_id']))
            else:
                #turn scientific_name from a list to a string
                ids_list = get_label(ids, site_key)
                datasets.append(ids_list)
                table = sets[1] #list of all ids in dataset
                samples.append(list(table['climb_id'])) #climb id

    # Using filter and lambda to remove 'mscape' from 'all' for all_other_directories list
    mscape_datasets = list(filter(lambda item: item not in public_datasets, datasets))

    # Read bacteria load
    data = make_microbial_count_table(reports, grouped_metadata, site_key)

    microbe_types = ["All", "Bacteria", "Fungi", "Viruses", "Archaea", "Protists"]
    load_list = []

    for microbe_type in microbe_types:
        data_sorted = data.sort_values(microbe_type)
        
        dataset = data_sorted['index']
        average = data_sorted[microbe_type]
        
        colors = []
        
        for name in dataset:
            # Find corresponding value in column B
            output = list(palette.loc[palette['name_order'] == name, 'colours'])
            output = ''.join(output)
            colors.append(output)

        # Figure Size
        fig, ax = plt.subplots(figsize =(8, 6))
        # convert y-axis to Logarithmic scale
        plt.xscale("log")
        # Horizontal Bar Plot
        ax.barh(dataset, average, color=colors)
        
        y = np.asarray([i for i in range(len(dataset))])
        ax.set_yticks(y)
        yid = np.asarray([i for i in dataset])
        ax.set_yticklabels(yid)
        
            
        new_labels = []
        
        for element in ax.yaxis.get_ticklabels():
            if "public" not in element.get_text().lower():
                element.set_fontweight('bold')
                new_labels.append(element)
            else:
                new_labels.append(element)
                
        ax.set_yticklabels(new_labels)
        plt.title(microbe_type)
        plt.xlabel('Average Read Count (Log)', fontweight='bold', horizontalalignment='center')
        
        
        # Save the figure as a Base64 string
        buf = io.BytesIO()
        plt.savefig(buf, format='png', bbox_inches='tight')
        buf.seek(0)
        microbe_load = base64.b64encode(buf.read()).decode('utf-8')
        load_list.append(microbe_load)
        buf.close()
        plt.close(fig)

    #shannon plots from R (in html code)
    with open(plots_dir + 'total_diversity.png', "rb") as image_file:
        total_diversity = base64.b64encode(image_file.read()).decode('utf-8')
        
    #shannon plots from R (in html code)
    with open(plots_dir + 'dna_diversity.png', "rb") as image_file:
        dna_diversity = base64.b64encode(image_file.read()).decode('utf-8')

    #shannon plots from R (in html code)
    with open(plots_dir + 'rna_diversity.png', "rb") as image_file:
        rna_diversity = base64.b64encode(image_file.read()).decode('utf-8')

    #shannon plots from R (in html code)
    with open(plots_dir + 'total_evenness.png', "rb") as image_file:
        total_evenness = base64.b64encode(image_file.read()).decode('utf-8')
        
    #shannon plots from R (in html code)
    with open(plots_dir + 'dna_evenness.png', "rb") as image_file:
        dna_evenness = base64.b64encode(image_file.read()).decode('utf-8')

    #shannon plots from R (in html code)
    with open(plots_dir + 'rna_evenness.png', "rb") as image_file:
        rna_evenness = base64.b64encode(image_file.read()).decode('utf-8')
        
    # get microbe taxa
    def make_two_sets(reports, grouped_metadata, column_to_drop, filter_count):
        taxons = ['F', 'G', 'S']
        richness_tables = []
        
        for taxon in taxons:
            df = make_richness_table(reports, grouped_metadata, taxon, filter_count, site_key)
            df = df.drop(columns=[column_to_drop])

            numerical_df = df.drop(columns=["index"])
            df["all"] = numerical_df.sum(axis=1)
            df_sorted = df.sort_values("all")
            df_sorted = df_sorted.drop(columns=["all"])
            
            if taxon == "F":
                word = "Families"
            elif taxon == "G":
                word = "Genera"
            else:
                word = "Species"
                
            # Generate the plot
            fig, ax = plt.subplots(figsize=(8, 6))
            #https://davidmathlogic.com/colorblind/#%237F9DEA-%23785EF0-%23DC267F-%23F58744-%23FFD900
            colors = ['#7F9CEA', '#785EF0', '#DC267F', '#F58744', '#FFD900']
            df_sorted.plot(
                x='index', kind='bar', stacked=True,
                title='Number of Microbial '+word, ax=ax,
                color=colors
            )

            dataset = df_sorted['index']
            x = np.asarray([i for i in range(len(dataset))])
            ax.set_xticks(x)
            xid = np.asarray([i for i in dataset])
            ax.set_xticklabels(xid)
            
                
            new_labels = []
            
            for element in ax.xaxis.get_ticklabels():
                if "public" not in element.get_text().lower():
                    element.set_fontweight('bold')
                    new_labels.append(element)
                else:
                    new_labels.append(element)
                    
            ax.set_xticklabels(new_labels)
            
            # Save the figure as a Base64 string
            buf = io.BytesIO()
            plt.savefig(buf, format='png', bbox_inches='tight')
            buf.seek(0)
            richness_table = base64.b64encode(buf.read()).decode('utf-8')
            richness_tables.append(richness_table)
            buf.close()
            plt.close(fig)

        return richness_tables

    filter_count = 5
    all_bac_tables = make_two_sets(reports, grouped_metadata, f"Bacteria > {filter_count}", filter_count)
    filtered_bac_tables = make_two_sets(reports, grouped_metadata, "Bacteria", filter_count)

    #make bar plots for total read counts (for heatmaps)
    sns.set_style("white")
    #heatmap function
    average, mscape, mscape_names, both_counts = get_heatmap(reports, grouped_metadata, site_key)

    mscape_colours = []
    for name in mscape_names:
        # Find corresponding value in column "colours"
        output = list(palette.loc[palette['name_order'] == name, 'colours'])
        output = ''.join(output)
        mscape_colours.append(output)

    two_total_counts = []
    total_samples = 0

    #Get number of total samples
    for df in both_counts[0]:
            samples = df.index
            total_samples = total_samples + len(samples)

    # Get list of relative width ratios for each subplot
    width_ratios = []
    for df in both_counts[0]:
        samples = df.index
        width = len(samples)/total_samples
        width_ratios.append(width)

        no_mscape = len(both_counts[0])

    # Figure Size
    fig, ax = plt.subplots(figsize=(18,6), ncols=no_mscape, gridspec_kw={'width_ratios': width_ratios})

    # Normalise counts for all mscape dataset subplots on the same y-scale
    all_counts = []

    for count_list in both_counts[0]: 
        counts = count_list.tolist() 
        all_counts = all_counts + counts

    # Convert strings to integers using
    # list comprehension
    int_counts = [int(item) for item in all_counts]

    max_value = max(int_counts)
    upper_lim = max_value + (max_value/10) # ensures linear graphs will be 110% in height/y-axis

    #make all mscape subplots
    loop = 0
    for count_list in both_counts[0]:
        df = pd.DataFrame({"Total Counts": count_list})
        df["Total Counts"] = pd.to_numeric(df["Total Counts"])    
        counts = df["Total Counts"]
        index = df.index.to_list()
        # Convert strings to integers using
        # list comprehension
        int_index = [int(item) for item in index]

        ax = plt.subplot(1, no_mscape, loop+1)

        # Optional: convert y-axis to Logarithmic scale
        #plt.yscale("log")
        # OR: set y-lim to something above 10% of highest value
        ax.set_ylim([0, upper_lim])

        # Horizontal Bar Plot
        ax.bar(int_index, counts, color=mscape_colours[loop])
        ax.ticklabel_format(style='plain')
        plt.xticks([])

        if loop == 0:
            plt.ylabel('Total Read Count', fontweight='bold', ha='center', labelpad=20)
            plt.xlabel(f'{mscape_names[loop]}', fontweight='bold', horizontalalignment='center')
        elif loop < (no_mscape-1):
            plt.yticks([])      
            plt.xlabel(f'{mscape_names[loop]}', fontweight='bold', horizontalalignment='center')
        elif loop == (no_mscape-1):
            plt.yticks([])
            plt.xlabel(f'{mscape_names[loop]}', fontweight='bold', horizontalalignment='center')
        loop += 1
        plt.subplots_adjust(wspace=0.02, hspace=0)
        plt.xlim([0,len(int_index)])

    # Save the figure as a Base64 string
    buf = io.BytesIO()
    plt.savefig(buf, format='png', bbox_inches='tight')
    buf.seek(0)
    counts = base64.b64encode(buf.read()).decode('utf-8')
    two_total_counts.append(counts)
    buf.close()
    plt.close(fig)

    #make a bar graph for public dataset dataframe
    pub_df = both_counts[1]
    df = pd.DataFrame({"Total Counts": pub_df})
    df["Total Counts"] = pd.to_numeric(df["Total Counts"])    
    counts = df["Total Counts"]
    index = df.index.to_list()

    all_counts = both_counts[1].tolist()
    int_counts = [int(item) for item in all_counts]

    max_value = max(int_counts)
    upper_lim = max_value + (max_value/10)

    # Figure Size
    fig, ax = plt.subplots(figsize=(18,6))

    # optional: convert y-axis to Logarithmic scale
    #plt.yscale("log")
    # OR: set y-lim to something above 10% of highest value
    ax.set_ylim([0, upper_lim])
    # Horizontal Bar Plot
    ax.bar(index, counts, color="lightgray")
    ax.ticklabel_format(style='plain')

    plt.xticks([])
    plt.xlabel("Public data", fontweight='bold', horizontalalignment='center')
    plt.ylabel('Total Read Count', fontweight='bold', ha='center', labelpad=20) 

    plt.xlim([0,len(index)])

    # Save the figure as a Base64 string
    buf = io.BytesIO()
    plt.savefig(buf, format='png', bbox_inches='tight')
    buf.seek(0)
    counts = base64.b64encode(buf.read()).decode('utf-8')
    two_total_counts.append(counts)
    buf.close()
    plt.close(fig)
    

    def get_two_maps(sort_way):
        heatmap_list = []
        total_samples = 0
        for df in sort_way[0]:
            df.set_index('Scientific_Name', inplace=True)
            samples = df.columns.tolist()
            total_samples = total_samples + len(samples)

        width_ratios = []
        for df in sort_way[0]:
            samples = df.columns.tolist()
            width = len(samples)/total_samples
            width_ratios.append(width)

        no_mscape = len(average[0])

        fig, ax = plt.subplots(figsize=(18,6), ncols=no_mscape, gridspec_kw={'width_ratios': width_ratios})

        loop = 0
        for df in sort_way[0]:
            df.replace(0, np.nan, inplace=True)
            
            genus = df.index.tolist()
            samples = df.columns.tolist()
            matrix = df.to_numpy()
            
            heatmap = np.reshape(matrix, (len(genus), len(samples)))

            ax = plt.subplot(1, no_mscape, loop+1)   

            plot = ax.pcolormesh(samples, genus, heatmap, vmin=0, vmax=100, cmap='magma')
            
            plt.xticks([])

            if loop == 0:
                plt.yticks(genus, rotation=0,fontsize='10')
                plt.xlabel(mscape_names[loop], fontweight='bold', horizontalalignment='center')
            elif loop < (no_mscape-1):
                plt.yticks([])      
                plt.xlabel(mscape_names[loop], fontweight='bold', horizontalalignment='center')
            elif loop == (no_mscape-1):
                plt.yticks([])
                fig.colorbar(plot, ax=ax, aspect=30)
                plt.xlabel(mscape_names[loop], fontweight='bold', horizontalalignment='center')
            loop += 1
        plt.subplots_adjust(wspace=0.02, hspace=0)

        # Save the figure as a Base64 string
        buf = io.BytesIO()
        plt.savefig(buf, format='png', bbox_inches='tight')
        buf.seek(0)
        heatmap = base64.b64encode(buf.read()).decode('utf-8')
        heatmap_list.append(heatmap)
        buf.close()
        plt.close(fig)
        
        pub_df = sort_way[1]
        pub_df.set_index('Scientific_Name', inplace=True)
        pub_df.replace(0, np.nan, inplace=True)
            
        genus = pub_df.index.tolist()
        samples = pub_df.columns.tolist()
        matrix = pub_df.to_numpy()
            
        heatmap = np.reshape(matrix, (len(genus), len(samples))) 

        fig, ax = plt.subplots(figsize=(18.5,5.5))
        plot = ax.pcolormesh(samples, genus, heatmap, vmin=0, vmax=100, cmap='magma')
        fig.colorbar(plot, ax=ax, aspect=30, pad=0.01)

        plt.xticks([])
        plt.yticks(genus, rotation=0,fontsize='10')
        plt.xlabel("Public data", fontweight='bold', horizontalalignment='center')
        plt.subplots_adjust(wspace=0.02, hspace=0) 
        
        # Save the figure as a Base64 string
        buf = io.BytesIO()
        plt.savefig(buf, format='png', bbox_inches='tight')
        buf.seek(0)
        heatmap = base64.b64encode(buf.read()).decode('utf-8')
        heatmap_list.append(heatmap)
        buf.close()
        plt.close(fig)

        return heatmap_list

    # Make two lists, sorted by mscape and by all
    mscape_maps = get_two_maps(mscape)
    average_maps = get_two_maps(average)

    output_path = args.final_report
    os.makedirs(output_path, exist_ok=True) #make directory path into real directory

    template_dir = args.template
    # Render the template with the Base64 string
    template = Template(filename=template_dir)
    html_content = template.render(all_load=load_list[0],
                                bac_load=load_list[1], 
                                fungi_load=load_list[2],
                                virus_load=load_list[3],
                                archaea_load=load_list[4],
                                protist_load=load_list[5],

                                filter_count=filter_count,
                                all_family_richness=all_bac_tables[0], all_genus_richness=all_bac_tables[1], all_species_richness=all_bac_tables[2],
                                f_family_richness=filtered_bac_tables[0], f_genus_richness=filtered_bac_tables[1], f_species_richness=filtered_bac_tables[2],
                                
                                total_diversity=total_diversity, total_evenness=total_evenness,
                                dna_diversity=dna_diversity, dna_evenness=dna_evenness,
                                rna_diversity=rna_diversity, rna_evenness=rna_evenness,
                                mscape_counts=two_total_counts[0], public_counts=two_total_counts[1],
                                mscape_by_total=average_maps[0], public_by_total=average_maps[1],
                                mscape_by_mscape=mscape_maps[0], public_by_mscape=mscape_maps[1])

    # Save the rendered HTML to a file
    with open(f"{output_path}negcontm_summary.html", "w") as f:
        f.write(html_content)

    print(f"HTML file generated: {output_path}")

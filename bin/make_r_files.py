#!/usr/bin/env python

import os
import pandas as pd
import numpy as np
import argparse

spikeins = {
        12242:["Tobamovirus","Tobacco_mosaic_virus"]
    }

def get_broad_count(needed_samples, reports, microbe_type, taxon_level):

    # Initialize an empty list to store dataframes
    dfs = []
    # Loop over each kraken file in reports directory
    for sample in needed_samples:
        for filename in os.listdir(reports):
            if f'{sample}.' in filename and filename.endswith('report.txt'):
                #open new file and read it line by line
                file = open(reports+"/"+filename)
        
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
                    if "D" in current_rank:
                        start_reading = False
                        if microbe_type == "All":
                            start_reading = True #read everything
                        elif microbe_type == "DNA":
                            if current_name != "Riboviria": #read everything except rna
                                start_reading = True
                        elif microbe_type == "RNA":
                            if current_name == "Riboviria": #only read rna
                                start_reading = True
                    elif start_reading:
                        read_counts.append(line.split()[1])
                        rank.append(line.split()[3])            
                        #this turns scientific name (which sometimes have multiple words) into a list within a list
                        sci_name.append(line.split()[5:])


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
    
    
    # Define a function to convert a column to numeric type if it's not 'rank' or 'scientific name'
    def convert_to_numeric(column):
        if column.name not in ['Rank', 'Scientific_Name']:
            return pd.to_numeric(column, errors='coerce')
        else:
            return column
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
    
    # Boolean indexing to filter rows showing only taxon wanted in the rank column
    filtered_df = numeric_df[numeric_df[columns_to_search].apply(lambda x: x.isin(keywords).any(), axis=1)]
    
    empty_data = {'Scientific_Name': 0}
    transposed_df = pd.DataFrame(empty_data, index=['Scientific_Name'])
    richness_df = pd.DataFrame([])
    
    if len(filtered_df.index) > 0:
        # Rearrange column 'Scientific Name' to the first position
        wordy_columns = ['Scientific_Name', 'Rank']
        # check if columns are not in the wordy_columns list
        column_order = ['Scientific_Name'] + ['Rank'] + [col for col in filtered_df.columns if col not in wordy_columns]
        filtered_df = filtered_df[column_order]
        
        # Turning it into percentages
        columns_to_sum = filtered_df.drop(columns=wordy_columns)   
        sample_sums = columns_to_sum.sum(axis=0)
        # Define a function to convert a column to % within sample if it's not 'rank' or 'scientific name'
        def percentify(column):
            if column.name not in ['Scientific_Name', 'Rank']:
                name = column.name
                return column/sample_sums[name]
            else:
                return column
        # Apply the function to each column
        percent_df = filtered_df.apply(percentify)
        drop_df = percent_df.drop(columns=["Rank"])
        drop_df.fillna(value=0.0, inplace=True)
        drop_df.set_index('Scientific_Name', inplace=True)
        transposed_df = drop_df.transpose()
        
        number = filtered_df.shape[0]
            
        data = {'Name': microbe_type, 'Number of Genera': number}
        richness_df = pd.DataFrame(data, index=[0])
        richness_df = richness_df.fillna(value=np.nan)
        
    return richness_df, transposed_df
        
# make 3 lists of richness/diversity by DNA, RNA, and all combined
def get_broad_category(set, needed_samples, reports, taxon_level):
    taxa = ["All", "DNA", "RNA"]
    current_dfs = []
    percent_dfs = []
    for microbe_type in taxa:
        rich, perc = get_broad_count(needed_samples, reports, microbe_type, taxon_level)
        current_dfs.append(rich)
        percent_dfs.append(perc)
        
    #merge all dataframes together for richness table
    table = pd.concat(current_dfs, axis = 0, join = "outer")
    table = table.groupby('Name', as_index=True).first()
    table.fillna(value=0.0, inplace=True)

    transposed_df = table.transpose()
    transposed_df.reset_index(inplace=True)

    transposed_df['index'] = set
    return transposed_df, percent_dfs

# Do the above for all directories
def make_richness_table(reports, grouped_metadata, taxon_level):
    #group by site
    datasets = []
    samples = []

    for sets in grouped_metadata:
        ids = list(sets[0])
        if "public" in ids[0].lower(): #if dataset is public
            subgroup = sets[1].groupby(['run_id'])
            for subset in subgroup:
                run_id = ''.join(subset[0])
                #turn scientific_name from a list to a string
                ids_list = '_'.join(ids)
                ids_list = ids_list.replace('public', 'public_'+run_id) #add run_id to name
                ids_list = ids_list.replace('_extraction_control', '')
                ids_list = ids_list.replace('_resp_matrix_mc110', '_matrix')

                datasets.append(ids_list)
                table = subset[1] #list of all ids in sub_dataset
                samples.append(list(table['climb_id'])) #climb id
        else:
            #turn scientific_name from a list to a string
            ids_list = '_'.join(ids)
            ids_list = ids_list.replace('_extraction_control', '')
            ids_list = ids_list.replace('_resp_matrix_mc110', '_matrix')
            datasets.append(ids_list)
            table = sets[1] #list of all ids in dataset
            samples.append(list(table['climb_id'])) #climb id

    loop = 0
    single_dfs = []
    single_perc_dfs = []
    for set in datasets:
        needed_samples = samples[loop]
        transposed_df, percent_dfs = get_broad_category(set, needed_samples, reports, taxon_level)
        single_dfs.append(transposed_df)
        single_perc_dfs.append(percent_dfs)
        loop += 1
        
    #merge all dataframes together
    final_table = pd.concat(single_dfs, axis = 0, join = "outer")
    final_table = final_table.groupby('index', as_index=False).first()
    final_table.fillna(value=0.0, inplace=True)

    return final_table, single_perc_dfs 

reports = "/Users/angelasun/Downloads/neg_control_reports/all_reports"
metadata = pd.read_csv("/Users/angelasun/Downloads/neg_control_reports/public_metadata.csv")
    
#make output text_files directory
output_dir = "/Users/angelasun/Downloads/local-test/test_files/"
os.makedirs(output_dir, exist_ok=True) #make directory path into real directory
#filter for richness data on genus level
taxon = "G"

grouped_metadata = metadata.groupby(['site', 'control_type_details'])
richness, diversity = make_richness_table(reports, grouped_metadata, taxon)

richness.to_csv(output_dir+"richness_table.txt", sep='\t', index=False)

names = []
for sets in grouped_metadata:
    ids = list(sets[0])
    if "public" in ids[0].lower(): #if dataset is public
        subgroup = sets[1].groupby(['run_id'])
        for subset in subgroup:
            run_id = ''.join(subset[0])
            #turn scientific_name from a list to a string
            ids_list = '_'.join(ids)
            ids_list = ids_list.replace('public', 'public_'+run_id) #add run_id to name
            ids_list = ids_list.replace('_extraction_control', '')
            ids_list = ids_list.replace('_resp_matrix_mc110', '_matrix')
            names.append(ids_list)
    else:
        #turn scientific_name from a list to a string
        ids_list = '_'.join(ids)
        ids_list = ids_list.replace('_extraction_control', '')
        ids_list = ids_list.replace('_resp_matrix_mc110', '_matrix')
        names.append(ids_list)

#save all text_files in working/processing output directory
loop = 0
for df_set in diversity:
    df_set[0].to_csv(os.path.join(output_dir, f"{names[loop]}.total.txt"), sep='\t', index=False)
    df_set[1].to_csv(os.path.join(output_dir, f"{names[loop]}.dna.txt"), sep='\t', index=False)
    df_set[2].to_csv(os.path.join(output_dir, f"{names[loop]}.rna.txt"), sep='\t', index=False)
    loop += 1


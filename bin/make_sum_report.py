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
from natsort import natsorted
import json

from utils import spikeins, get_label, make_count_and_perc_dfs, define_datasets

#Merge all taxon count dataframes together to get all microbe type count data per dataset
def get_microbial_load(set, needed_samples, reports):
    count_df = make_count_and_perc_dfs(needed_samples, reports, "count")
    #Get list of spikeins and make a spikeins dataframe from main dataframe
    spikein_names = []
    for spike in spikeins:
        spike_list = spikeins[spike]
        spikein_names.append(spike_list[1])

    if isinstance(spikein_names, str):
        spike_df = count_df.loc[count_df["Scientific_Name"] == spikein_names]
    else:
        spike_df = count_df.loc[count_df["Scientific_Name"].isin(spikein_names)]

    #Get sum of spikeins
    spike_df = spike_df.drop(columns=['Scientific_Name', 'Rank', "Domain"])
    copy_df = spike_df.copy()
    copy_df["Sum"] = copy_df.sum(axis=1) #sum for each spikein taxon
    spike_sum = copy_df["Sum"].sum(axis=0) #sum for all spikein taxa combined

    #Find total counts for microbe categories that can be found in single rows
    domains_list = ["Bacteria", "Archaea", "Fungi", "Viruses"]
    domains_df = count_df.loc[count_df["Scientific_Name"].isin(domains_list)]
    domains_numeric = domains_df.drop(columns=["Scientific_Name", "Rank", "Domain"])
    domains_df["Sum"] = domains_numeric.sum(axis=1)

    #Make initial table
    table_columns = ["Scientific_Name", "Sum"]
    current_df = domains_df[table_columns].copy()
    current_df = current_df.reset_index(drop=True)

    #Remove spikein from virus sum
    virus_index = current_df.index.get_loc(current_df[current_df["Scientific_Name"] == "Viruses"].index[0])
    virus_sum = current_df.loc[virus_index, "Sum"]
    new_virus_sum = virus_sum - spike_sum
    current_df.loc[virus_index, "Sum"] = new_virus_sum

    #Find Protists sum by making a protist df
    protist_list = ["Sar", "Discoba", "Metamonada"]
    protist_df = count_df.loc[count_df["Scientific_Name"].isin(protist_list)]
    protist_numeric = protist_df.drop(columns=["Scientific_Name", "Rank", "Domain"])
    protist_df["Sum"] = protist_numeric.sum(axis=1)
    total_protists = protist_df["Sum"].sum(axis=0)

    # Create a dictionary with the data for the new row
    protist_row = {'Scientific_Name': 'Protists', 'Sum': total_protists}
    # Append the dictionary to the DataFrame
    current_df.loc[len(current_df)] = protist_row

    #Find total sum
    total_sum = current_df["Sum"].sum()
    all_row = {'Scientific_Name': 'All', 'Sum': total_sum}
    # Append the dictionary to the DataFrame
    current_df.loc[len(current_df)] = all_row

    current_df.set_index("Scientific_Name", inplace=True)

    transposed_df = current_df.transpose()
    transposed_df.reset_index(inplace=True)

    transposed_df['index'] = set
    return transposed_df

#Make microbe count table for each dataset and merge them together
def make_microbial_count_table(reports, grouped_metadata, site_key):

    datasets, samples = define_datasets(grouped_metadata, site_key)

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

def get_all_taxa(set, needed_samples, reports, taxon_level, filter_count):
    count_df = make_count_and_perc_dfs(needed_samples, reports, "count")
    
    #remove spikeins
    for spike in spikeins:
        count_df = count_df.loc[~count_df["Scientific_Name"].astype(str).isin(spikeins[spike])]

    # Define keywords and columns to search
    keywords = [taxon_level] #the taxonomy level I want to filter for
    columns_to_search = ['Rank']
    
    # Boolean indexing to filter rows showing only genus in the rank column
    filtered_df = count_df[count_df[columns_to_search].apply(lambda x: x.isin(keywords).any(), axis=1)]

    if len(filtered_df.index) > 0 :
        taxa = ["Bacteria", f"Bacteria > {filter_count}", "Fungi", "Viruses", "Archaea", "Protists"]
        taxon_counts = []
        for taxon in taxa:
            if taxon == f"Bacteria > {filter_count}":
                current_taxa_df = filtered_df.loc[filtered_df["Domain"] == "Bacteria"]
                numeric_columns = current_taxa_df.drop(columns=['Scientific_Name', 'Rank', 'Domain'])
                column_sums = numeric_columns.sum(axis=1) #sum of counts
                no_columns = numeric_columns.shape[1] #number of columns
                current_taxa_df["Average"] = column_sums/no_columns
                filtered_count = current_taxa_df[current_taxa_df["Average"] >= filter_count]
                taxon_count = filtered_count.shape[0]

            else:
                current_taxa_df = filtered_df.loc[filtered_df["Domain"] == taxon]
                taxon_count = current_taxa_df.shape[0]
            
            taxon_counts.append(taxon_count)

    new_df = pd.DataFrame({"Name": taxa, set: taxon_counts})
    new_df.set_index("Name", inplace=True)
    transposed_df = new_df.transpose()
    transposed_df.reset_index(inplace=True)
    return transposed_df

def make_richness_table(reports, grouped_metadata, taxon_level, filter_count, site_key):
 
    datasets, samples = define_datasets(grouped_metadata, site_key)

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
    perc_df = make_count_and_perc_dfs(needed_samples, reports, "perc")
   
    # Select columns to sum (excluding 'Rank' and 'Scientific Name')
    columns_to_sum = perc_df.drop(columns=['Rank', 'Scientific_Name', 'Domain'])
    # Calculate the sum of values in each row
    column_sums = columns_to_sum.sum(axis=1)
    # Count the number of columns
    num_columns = perc_df.shape[1]
    #either get a sum of all seqs, or normalise it by dividing by num_columns
    perc_df["Perc_Seqs_Overall"] = column_sums/num_columns
    
    #remove spikeins
    for spike in spikeins:
        perc_df = perc_df.loc[~perc_df["Scientific_Name"].astype(str).isin(spikeins[spike])]

    #Filtering the dataframe    
    # Define keywords and columns to search
    keywords = ['G'] #the taxonomy level I want to filter for
    columns_to_search = ['Rank']
    
    # Boolean indexing to filter rows showing only genus in the rank column
    genus_df = perc_df[perc_df[columns_to_search].apply(lambda x: x.isin(keywords).any(), axis=1)]
    
    # and remove rows with zero values (0 or 0.0)
    keywords_to_remove = [0]
    # Boolean indexing to filter rows in % of seqs overall that contain 0
    filtered_df = genus_df[~genus_df["Perc_Seqs_Overall"].astype(float).isin(keywords_to_remove)]
    
    # Rearrange column 'Scientific Name' to the first position
    wordy_columns = ['Scientific_Name']
    # check if columns are not in the wordy_columns list
    column_order = ['Scientific_Name'] + [col for col in filtered_df.columns if col not in wordy_columns]
    filtered_df = filtered_df[column_order]
    
    #drop column “rank” since they are all Gs and drop "domain"
    no_rank_df = filtered_df.drop(columns=['Rank', 'Domain'])
    return no_rank_df

def get_heatmap(reports, grouped_metadata, site_key):
    
    #group by site
    datasets = []
    public_datasets = []
    samples = []
    public_samples = []
    sample_dates = []
    all_samples = ["Domain", "Scientific_Name", "Rank"]

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
                    dates_columns = ['collection_date', 'received_date']
                    sample_dates.append(table[dates_columns]) #sample dates

                    all_samples.extend(list(table['climb_id']))
            else:
                #turn scientific_name from a list to a string
                ids_list = get_label(ids, site_key)
                datasets.append(ids_list)
                table = sets[1] #list of all ids in dataset
                samples.append(list(table['climb_id'])) #climb id
                dates_columns = ['collection_date', 'received_date']
                sample_dates.append(table[dates_columns]) #sample dates

                all_samples.extend(list(table['climb_id']))
            
    
    average_list = []
    og_list = []

    loop = 0
    for set in datasets:
        needed_samples = samples[loop] #our current set of sample names
        perc_df, og_df = make_perc_df(needed_samples, reports)

        sample_date_df = sample_dates[loop]

        date_list = ['not_applicable']
        for index, row in sample_date_df.iterrows():
            if row['collection_date'] == row['collection_date']:
                date_list.append(row['collection_date'])
            elif row['received_date'] == row['received_date']:
                date_list.append(row['received_date'])
            else:
                date_list.append('NaN')
        date_list.append(100)

        perc_df = perc_df.rename(columns={c: f"{set}_"+c for c in perc_df.columns if c not in ['Scientific_Name', 'Perc_Seqs_Overall']})
        perc_df.loc[len(perc_df)] = date_list

        row_number = perc_df.index.get_loc(perc_df[perc_df["Scientific_Name"] == "not_applicable"].index[0])
        perc_df.iloc[row_number, 1:-1] = pd.to_datetime(perc_df.iloc[row_number, 1:-1], format='%d/%m/%Y')
        
        #Change all "average percentage" columns to their respective dataset names to avoid clashes when merging
        perc_df[set] = perc_df["Perc_Seqs_Overall"]
        perc_df = perc_df.drop(columns=["Perc_Seqs_Overall"])
        average_list.append(perc_df)

        og_list.append(og_df)
        loop += 1

    #First define the merged dataframe by the starting dataframe
    loop_count = 0
    merge_df = average_list[0]
    og_merge_df = og_list[0]
    #Then merge all other dataframes to the pre-existing dataframe
    for df in average_list:
        loop_count += 1
        if loop_count < len(datasets):
            current_df = average_list[loop_count]
            merge_df = merge_df.merge(current_df, on="Scientific_Name", how="outer")
            current_og_df = og_list[loop_count]
            og_merge_df = og_merge_df.merge(current_og_df, on=["Scientific_Name", "Rank", "Domain"], how="outer")
    
    merge_df.fillna(value=0, inplace=True)
    og_merge_df.fillna(value=0, inplace=True)

    # Rearrange column 'Scientific Name' to the first position
    wordy_columns = ['Scientific_Name', 'Rank', 'Domain']
    # check if columns are not in the wordy_columns list
    column_order = ['Domain'] + ['Scientific_Name'] + ['Rank'] + [col for col in og_merge_df.columns if col not in wordy_columns]
    og_merge_df = og_merge_df[column_order]
    og_merge_df.columns = all_samples
    
    os.makedirs(f'{output_path}dataframes/', exist_ok=True)
    og_merge_df.to_csv(f'{output_path}dataframes/percentage_df.csv', index=False)

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
        both_dates = []

        mscape_dfs = []
        publics = "public"

        #make specific dfs for all dataframes in mscape_datasets
        for dataset_name in mscape_datasets:
            short_name = dataset_name.split("(")[0]
            category = dataset_name.split("(")[1].replace(")","")
            name_col = "Scientific_Name"
            df = total_df.loc[:, total_df.columns.str.contains(short_name) & total_df.columns.str.contains(category)]
            df[name_col] = total_df[name_col]
            mscape_dfs.append(df)
        
        loop_count = 0
        sorted_mscapes = []
        sorted_m_counts = []
        sorted_m_dates = []
        for mscape_df in mscape_dfs:
            #Get rid of dataset as prefixes
            mscape_df = mscape_df.rename(columns={c: c.replace(f"{mscape_datasets[loop_count]}_", "") for c in mscape_df.columns if c not in ['Scientific_Name']})
            row_number = mscape_df.index.get_loc(mscape_df[mscape_df["Scientific_Name"] == "not_applicable"].index[0])
            mscape_df.iloc[row_number] = natsorted(mscape_df.iloc[row_number])
            new_row_number = mscape_df.index.get_loc(mscape_df[mscape_df["Scientific_Name"] == "not_applicable"].index[0])
            
            dates = list(mscape_df.iloc[new_row_number])
            sorted_m_dates.append(dates)
            mscape_df = mscape_df[~mscape_df.Scientific_Name.str.contains("not_applicable")]          
            sorted_mscapes.append(mscape_df)

            mscape_matrix = mscape_df.drop(columns=["Scientific_Name"])
            mscape_counts = mscape_matrix.transpose()
            mscape_counts = mscape_counts.reset_index()
            sorted_m_counts.append(mscape_counts['index'])
            loop_count += 1
        
        public_df = total_df.loc[:, total_df.columns.str.contains(publics, case=False)]
        public_df["Scientific_Name"] = total_df["Scientific_Name"]
        public_df = public_df.rename(columns={c: c.split("_")[-1] for c in public_df.columns if c not in ['Scientific_Name']})
        row_number = public_df.index.get_loc(public_df[public_df["Scientific_Name"] == "not_applicable"].index[0])
        public_df.iloc[row_number] = natsorted(public_df.iloc[row_number])

        public_dates = list(public_df.iloc[row_number])
        public_df = public_df[~public_df.Scientific_Name.str.contains("not_applicable")]

        public_matrix = public_df.drop(columns=["Scientific_Name"])
        public_counts = public_matrix.transpose()
        public_counts = public_counts.reset_index()

        both_dfs.append(sorted_mscapes)
        both_dfs.append(public_df)

        both_counts.append(sorted_m_counts)
        both_counts.append(public_counts['index'])

        both_dates.append(sorted_m_dates)
        both_dates.append(public_dates)

        return both_dfs, both_counts, both_dates

    sort_by_average, both_counts, both_dates = sorting(merge_df, samples, datasets, datasets)
    sort_by_mscape, both_counts, both_dates = sorting(merge_df, mscape_samples, mscape_datasets, datasets)

    return sort_by_average, sort_by_mscape, mscape_datasets, both_counts, both_dates

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse text files and a CSV file.")
    parser.add_argument('--reports', nargs='+', help="List of text files", required=True)
    parser.add_argument('--metadata', help="CSV file path", required=True)
    parser.add_argument('--site_key', help="JSON file specifying site name to number", required=True)
    parser.add_argument('--plots_dir', help="Shannon plots directory", required=True)
    parser.add_argument('--final_reports', help="Output directory", required=True)
    parser.add_argument('--template', help="HTMl template", required=True)
    args = parser.parse_args()

    reports = args.reports
    metadata = pd.read_csv(args.metadata)
    plots_dir = args.plots_dir # Import R plots
    palette = pd.read_csv(f'{plots_dir}/colour_palette.txt', delimiter='\t')

    output_path = args.final_reports
    os.makedirs(output_path, exist_ok=True) #make directory path into real directory

    sns.set_style("whitegrid")
    all_directories = []
    mscape_directories = []
    #group by site and other identifiers
    grouped_metadata = metadata.groupby(['site', 'control_type_details'])
    #site name to number
    site_key = {}
    with open(args.site_key, 'r') as f:
        site_key = json.load(f)

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
    average, mscape, mscape_names, both_counts, both_dates = get_heatmap(reports, grouped_metadata, site_key)

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
    
    #make heatmaps
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

    os.makedirs(f'{output_path}/summary_report/', exist_ok=True)
    # Save the rendered HTML to a file
    with open(f"{output_path}summary_report/neg_contamn_report.html", "w") as f:
        f.write(html_content)

    print(f"HTML file generated: {output_path}")

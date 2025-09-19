#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
import io
import base64
from mako.template import Template
import numpy as np
import seaborn as sns
import os
import argparse
from natsort import natsorted
import json

from utils import spikeins, convert_to_numeric, make_count_and_perc_dfs, define_datasets, define_heatmap_datasets, make_heatmap_df, make_date_list, split_dfs, save_labelled_df

#Merge all taxon count dataframes together to get all microbe type count data per dataset
def get_microbial_load(dataset, needed_samples, reports):
    count_df = make_count_and_perc_dfs(needed_samples, reports, "count")
    count_df = count_df.drop(columns=["Taxon_ID"])

    #Get list of spikeins and make a spikeins dataframe from main dataframe
    spikein_names = []
    for spike in spikeins:
        spike_list = spikeins[spike]
        spikein_names.append(spike_list[2])

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
    if (current_df[current_df["Scientific_Name"] == "Viruses"].size > 0):
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

    transposed_df['index'] = dataset
    return transposed_df

#Make microbe count table for each dataset and merge them together
def make_microbial_count_table(reports, grouped_metadata, site_key):

    datasets, samples = define_datasets(grouped_metadata, site_key)

    loop = 0
    single_dfs = []
    for dataset in datasets:
        needed_samples = samples[loop] #our current set of sample names
        transposed_df = get_microbial_load(dataset, needed_samples, reports) #get dataset name, sample names, and all reports
        single_dfs.append(transposed_df)
        loop += 1
        
    #merge all dataframes together
    final_table = pd.concat(single_dfs, axis = 0, join = "outer")
    final_table = final_table.groupby('index', as_index=False).first()
    final_table.fillna(value=0.0, inplace=True)

    os.makedirs(f'{output_path}dataframes/', exist_ok=True)
    final_table.to_csv(f'{output_path}/dataframes/total_microbe_counts.csv', index=False)
    return final_table

def get_all_taxa(dataset, needed_samples, reports, taxon_level, filter_count):
    og_count_df = make_count_and_perc_dfs(needed_samples, reports, "count")
    count_df = og_count_df.drop(columns=["Taxon_ID"])

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

    new_df = pd.DataFrame({"Name": taxa, dataset: taxon_counts})
    new_df.set_index("Name", inplace=True)
    transposed_df = new_df.transpose()
    transposed_df.reset_index(inplace=True)
    return transposed_df, og_count_df

def make_richness_table(reports, grouped_metadata, taxon_level, filter_count, site_key):
 
    datasets, samples = define_datasets(grouped_metadata, site_key)

    loop = 0
    single_dfs = []
    og_count_dfs = []
    for dataset in datasets:
        needed_samples = samples[loop] #our current set of sample names
        transposed_df, og_count_df = get_all_taxa(dataset, needed_samples, reports, taxon_level, filter_count)
        single_dfs.append(transposed_df)
        og_count_dfs.append(og_count_df)
        loop += 1
    
    #First define the merged dataframe by the starting dataframe
    loop_count = 0
    merge_df = og_count_dfs[0]
    #Then merge all other dataframes to the pre-existing dataframe
    for df in og_count_dfs:
        loop_count += 1
        if loop_count < len(og_count_dfs):
            current_df = og_count_dfs[loop_count]
            merge_df = merge_df.merge(current_df, on=["Scientific_Name", "Rank", "Taxon_ID", "Domain"], how="outer")
    
    # Rearrange column 'Scientific Name' to the first position
    wordy_columns = ['Scientific_Name', 'Rank', 'Domain', 'Taxon_ID']
    # check if columns are not in the wordy_columns list
    column_order = ['Domain'] + ['Scientific_Name'] + ["Taxon_ID"] + ['Rank'] + [col for col in merge_df.columns if col not in wordy_columns]
    merge_df = merge_df[column_order]

    merge_df.fillna(value=0, inplace=True)
    merge_df.to_csv(f'{output_path}/dataframes/count_df.csv', index=False)
        
    #merge all dataframes together
    final_table = pd.concat(single_dfs, axis = 0, join = "outer")
    final_table = final_table.groupby('index', as_index=False).first()
    final_table.fillna(value=0.0, inplace=True)

    final_table.to_csv(f'{output_path}/dataframes/{taxon_level}_level_richness.csv', index=False)
    return final_table 


def get_heatmap(reports, grouped_metadata, site_key, condition, hcids):
    
    mscape_datasets, mscape_samples, mscape_dates = define_heatmap_datasets(grouped_metadata, site_key)

    all_samples = ["Domain", "Scientific_Name", "Taxon_ID", "Rank"]
    for sample_group in mscape_samples:
        all_samples = all_samples + sample_group         
    
    average_list = []
    og_list = []

    loop = 0
    for dataset in mscape_datasets:
        needed_samples = mscape_samples[loop] #our current set of sample names
        perc_df, og_df = make_heatmap_df(needed_samples, reports, condition)

        sample_date_df = mscape_dates[loop]
        date_list = make_date_list(sample_date_df)

        if condition == "perc":
            average = "Perc_Seqs_Overall"
        else:
            average = "Counts_Overall"
        perc_df = perc_df.rename(columns={c: f"{dataset}_"+c for c in perc_df.columns if c not in ['Domain', 'Scientific_Name', average]})
        perc_df.loc[len(perc_df)] = date_list

        row_number = perc_df.index.get_loc(perc_df[perc_df["Scientific_Name"] == "Date"].index[0])
        perc_df.iloc[row_number, 2:-1] = pd.to_datetime(perc_df.iloc[row_number, 2:-1], format='%Y-%m-%d')
        #Sort all samples by date
        perc_df.iloc[row_number, 2:-1] = natsorted(perc_df.iloc[row_number, 2:-1])

        #Change all "average percentage" columns to their respective dataset names to avoid clashes when merging
        perc_df[dataset] = perc_df[average]
        perc_df = perc_df.drop(columns=[average])
        average_list.append(perc_df)

        og_list.append(og_df)
        loop += 1


    # Merge the DataFrames on a specific column
    merge_df = pd.concat(average_list, axis = 0, join= "outer")  # Change join to 'outer' for outer join
    #merge on scientific name - this means that no scientific name is repeated if it's present in different dataframes
    merge_df = merge_df.groupby(['Domain', 'Scientific_Name'], as_index=False).first()
    merge_df.fillna(0, inplace=True)

    if condition == "perc":
        og_merge_df = pd.concat(og_list, axis = 0, join= "outer")  
        og_merge_df = og_merge_df.groupby(['Domain', 'Scientific_Name', "Taxon_ID", "Rank"], as_index=False).first()
        og_merge_df.fillna(0, inplace=True)
        og_merge_df.columns = all_samples #return column headings to purely sample names

        og_merge_df.to_csv(f'{output_path}/dataframes/percentage_df.csv', index=False)

    #get list of dates
    timestamp_df = merge_df.drop(columns=mscape_datasets)
    row_number = timestamp_df.index.get_loc(timestamp_df[timestamp_df["Scientific_Name"] == "Date"].index[0])
    date_df = timestamp_df.iloc[row_number]
   
    all_dates = list(date_df)
    merge_df = merge_df[~merge_df.Scientific_Name.str.contains("Date")]

    if condition == "perc":
        # Select columns to sum ('Scientific Name')
        columns_to_sum = merge_df[mscape_datasets]

        # Calculate the sum of values in each row
        column_sums = columns_to_sum.sum(axis=1)
        merge_df['Average'] = column_sums/len(mscape_datasets)

        # Calculate the max of values in each row
        # merge_df['Max'] = columns_to_sum.max(axis=1)

        sorted_df = merge_df.sort_values(by="Average", ascending=False)
        top_df = sorted_df.head(20)

        total_df = top_df.drop(columns=mscape_datasets)
        total_df = total_df.drop(columns=["Average"])

        total_df = total_df.sort_values(by="Domain", ascending=True)

        #save heatmap df
        save_labelled_df(total_df, all_dates, output_path, condition)

        #Split dataframes and make a subdataframe per dataset
        sorted_mscapes = split_dfs(mscape_datasets, total_df)   
        return sorted_mscapes, mscape_datasets
    
    elif condition == 'thresh':
        classified_count = merge_df.loc[merge_df["Scientific_Name"].isin(["root", "unclassified"])]
        classified_count = classified_count.drop(columns=mscape_datasets)

        # Rearrange column 'Scientific Name' to the first position
        wordy_columns = ['Scientific_Name', 'Domain']
        wordy_columns.extend(mscape_datasets)

        # check if columns are not in the wordy_columns list or "counts_overall" (denoted by mscape_datasets)
        numeric_columns = [col for col in merge_df.columns if col not in wordy_columns]
        
        bacteria_df = merge_df.loc[merge_df["Domain"] == "Bacteria"]
        bacteria_df = bacteria_df[(bacteria_df[numeric_columns] >= 500).any(axis=1)]

        other_domains = ["Archaea", "Eukaryota", "Viruses", "Fungi", "Protists"]
        other_df = merge_df.loc[merge_df["Domain"].isin(other_domains)]
        other_df = other_df[(other_df[numeric_columns] >= 50).any(axis=1)]

        thresh_df = pd.concat([bacteria_df,other_df], axis=0, join="outer")
        thresh_df = thresh_df.groupby(['Domain', 'Scientific_Name'], as_index=False).first()
        
        # Select columns to sum ('Scientific Name')
        columns_to_sum = thresh_df[mscape_datasets]
        # Calculate the sum of values in each row
        column_sums = columns_to_sum.sum(axis=1)
        thresh_df['Count_Sum'] = column_sums

        # Calculate the max of values in each row
        # merge_df['Max'] = columns_to_sum.max(axis=1)

        sorted_df = thresh_df.sort_values(by="Count_Sum", ascending=False)
        total_df = sorted_df.drop(columns=mscape_datasets)
        total_df = total_df.sort_values(by="Domain", ascending=True)
        total_df = total_df.drop(columns=["Count_Sum"])

        save_labelled_df(total_df, all_dates, output_path, "count")

        #Change name of sample columns to sample total count
        #total_df = total_df.rename(columns={c: c.split("[")[0].replace("]", "") for c in total_df.columns if c not in ['Scientific_Name', 'Domain']})
    
        top_df = sorted_df.head(20)
        top_df = top_df.sort_values(by="Domain", ascending=True)
        top_df = top_df.drop(columns=["Count_Sum"])
        top_df = top_df.drop(columns=mscape_datasets)
        #top_df = top_df.rename(columns={c: c.split("[")[0].replace("]", "") for c in top_df.columns if c not in ['Scientific_Name', 'Domain']})
    
        #Split dataframes and make a subdataframe per dataset
        total_mscapes = split_dfs(mscape_datasets, total_df)
        top_mscapes = split_dfs(mscape_datasets, top_df)
        class_mscapes = split_dfs(mscape_datasets, classified_count)

        return top_mscapes, total_mscapes, class_mscapes
    
    elif condition == "count":
        #keep count_df column order for listing samples
        hcid_df = pd.DataFrame(columns=merge_df.columns) 
        hcid_df = hcid_df.drop(columns=mscape_datasets)
        #find all taxa with hcid counts 
        tables = []
        hcid_min = []
        for table in hcids:
            sample_id = table.split("_hcid")[0] #unique id
            hcid_table = pd.read_csv(table)
            
            for site_id_name in hcid_df.columns:
                if site_id_name not in ["Domain", "Scientific_Name"]:

                    site_name = site_id_name.split("_")[0]
                    id_name_list = site_id_name.split("_")[1:]
                    id_name = '_'.join(id_name_list)

                    if sample_id == id_name:
                        hcid_min.append(hcid_table["min_count"])
                        hcid_table[f'{site_name}_{sample_id}'] = hcid_table["mapped_count"]
                        needed_columns = ["name", f'{site_name}_{sample_id}']
                        sub_table = hcid_table[needed_columns]
                       
                        tables.append(sub_table)
        
        if len(tables) == 0:
            hcid_df["Scientific_Name"] = hcid_table["name"]
        elif len(tables) > 0:
            if len(tables) == 1:
                merged_table = tables[0]
            else:
                # Merge the DataFrames on a specific column
                merged_table = pd.concat(tables, axis = 0, join= "outer")  # Change join to 'outer' for outer join
                #merge on hcid name - this means that no scientific name is repeated if it's present in different dataframes
                merged_table = merged_table.groupby('name', as_index=False).first()

            #move taxa names over to hcid_df
            hcid_df["Scientific_Name"] = merged_table["name"]

            #replace hcid_df blank columns with merged_table columns if there are counts 
            for sample in hcid_df.columns:
                if sample in list(merged_table.columns):
                    hcid_df[sample] = merged_table[sample]

        hcid_df.fillna(0, inplace=True)
        hcid_mscapes = split_dfs(mscape_datasets, hcid_df)
        
        return hcid_mscapes





if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse text files and a CSV file.")
    parser.add_argument('--reports', nargs='+', help="List of text files", required=True)
    parser.add_argument('--metadata', help="CSV file path", required=True)
    parser.add_argument('--site_key', help="JSON file specifying site name to number", required=True)
    parser.add_argument('--hcids', nargs='+', help="List of hcid csv files", required=True)
    parser.add_argument('--plots_dir', help="Shannon plots directory", required=True)
    parser.add_argument('--final_reports', help="Output directory", required=True)
    parser.add_argument('--template', help="HTMl template", required=True)
    args = parser.parse_args()

    reports = args.reports
    metadata = pd.read_csv(args.metadata)
    plots_dir = args.plots_dir # Import R plots
    palette = pd.read_csv(f'{plots_dir}/colour_palette.txt', delimiter='\t')
    hcids = args.hcids
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
    top_mscape_perc, mscape_names = get_heatmap(reports, grouped_metadata, site_key, "perc", hcids)
    top_mscape_count, total_mscape_count, class_counts = get_heatmap(reports, grouped_metadata, site_key, "thresh", hcids)
    hcid_count = get_heatmap(reports, grouped_metadata, site_key, "count", hcids)

    stacked_class = []
    
    def make_stacked_class(class_counts, mscape_names, bar_type):
        mscape_colours = []
        for name in mscape_names:
            # Find corresponding value in column "colours"
            output = list(palette.loc[palette['name_order'] == name, 'colours'])
            output = ''.join(output)
            mscape_colours.append(output)

        total_samples = 0

        #Get number of total samples
        for df in class_counts:
                columns = int(df.shape[1])
                samples = columns-2
                total_samples = total_samples + samples

        # Get list of relative width ratios for each subplot
        width_ratios = []
        for df in class_counts:
            samples = df.shape[1] - 2
            width = samples/total_samples
            width_ratios.append(width)

        no_mscape = len(class_counts)

        # Figure Size
        fig, ax = plt.subplots(figsize=(18,6), ncols=no_mscape, gridspec_kw={'width_ratios': width_ratios})

        #make all mscape subplots
        loop = 0
        for df in class_counts:
            df = df.drop(columns=["Domain"])
            df.set_index("Scientific_Name", inplace = True)
            tdf = df.transpose()
            tdf.reset_index(inplace=True)
            
            if bar_type == "relative":
                # Turning it into percentages
                columns_to_sum = tdf.drop(columns=["index"])   
                sample_sums = columns_to_sum.sum(axis=1)
                sample_sums.replace(0, np.nan, inplace=True)
                sample_sums.fillna(0.000000001)
                # Define a function to convert a column to % within sample if it's not 'index'
                def make_relative(column):
                    if column.name not in ['index']:
                        return column/sample_sums
                    else:
                        return column
                # Apply the function to each column
                relative_df = tdf.apply(make_relative)

            ax = plt.subplot(1, no_mscape, loop+1)

            # Convert strings to integers using
            # list comprehension
            #int_counts = [int(item) for item in all_counts]

            #max_value = max(int_counts)
            # upper_lim = max_value + (max_value/10) # ensures linear graphs will be 110% in height/y-axis
            colours = [mscape_colours[loop]]
            colours.append("#adb5bd")

            #make all mscape subplots
            if bar_type == "absolute":
                tdf.plot(
                        x='index', kind='bar', stacked=True, ax=ax, width=0.98, color=colours
                    )
            else:
                relative_df.plot(
                        x='index', kind='bar', stacked=True, ax=ax, width=0.98, color=colours
                    )
            #ax.ticklabel_format(style='plain')
            plt.xticks([])

            if loop == 0:
                if bar_type == "absolute":
                    plt.ylabel('Total Read Count', fontweight='bold', ha='center', labelpad=20)
                else:
                    plt.ylabel('Relative Read Count', fontweight='bold', ha='center', labelpad=20)
                plt.xlabel("")
                ax.get_legend().remove()
            elif loop < (no_mscape-1):
                plt.yticks([])
                plt.xlabel("")     
                ax.get_legend().remove()
            elif loop == (no_mscape-1):
                plt.yticks([])
                plt.xlabel("")
                ax.get_legend().remove()

                legend_elements = []
                name_loop = 0
                for color in mscape_colours:
                    legend_elements.append(Patch(facecolor=color, label=mscape_names[name_loop]))
                    name_loop += 1

                legend_elements.append(Patch(facecolor="#adb5bd", label="unclassified"))
                # Create the figure
                ax.legend(handles=legend_elements, bbox_to_anchor=(1.04, 1), loc="upper left")
            loop += 1
            plt.subplots_adjust(wspace=0.02, hspace=0)
            #plt.xlim([0,len(int_index)])

        # Save the figure as a Base64 string
        buf = io.BytesIO()
        plt.savefig(buf, format='png', bbox_inches='tight')
        buf.seek(0)
        counts = base64.b64encode(buf.read()).decode('utf-8')
        stacked_class.append(counts)
        buf.close()
        plt.close(fig)
    
    make_stacked_class(class_counts, mscape_names, "absolute")
    make_stacked_class(class_counts, mscape_names, "relative")

    sns.set_style("dark")
    
    hcid_display = "display: none"
    heatmap_list = []
    hcid_list = []
    #make heatmaps
    def plot_heatmaps(dataframes, map_type):
        total_samples = 0
        for df in dataframes:
            df.set_index('Scientific_Name', inplace=True)
            if "Domain" in df.columns:
                df = df.drop(columns=["Domain"])
            samples = df.columns.tolist()
            total_samples = total_samples + len(samples)

        width_ratios = []
        for df in dataframes:
            samples = df.columns.tolist()
            width = len(samples)/total_samples
            width_ratios.append(width)

        df_max = []
        df_min = []
        
        for df in dataframes:
            df = df.drop(columns=["Domain"])
            df_max.append(max(df.max()))
            df_min.append(min(df.min()))

        no_mscape = len(dataframes)

        if map_type == "perc" or map_type == "top count":
            fig, ax = plt.subplots(figsize=(18,6), ncols=no_mscape, gridspec_kw={'width_ratios': width_ratios})
        else:
            height = len(dataframes[0].index)
            height = height * 0.3
            fig, ax = plt.subplots(figsize=(18, height), ncols=no_mscape, gridspec_kw={'width_ratios': width_ratios})


        loop = 0
        for df in dataframes:
            df.replace(0, np.nan, inplace=True)
            domain_list = list(df["Domain"])
            df = df.drop(columns=["Domain"])
            

            taxa = df.index.tolist()
            samples = df.columns.tolist()
            matrix = df.to_numpy()

            heatmap = np.reshape(matrix, (len(taxa), len(samples)))
            heatmap = np.array(heatmap, dtype=np.float64)  # Ensure it's numeric

            ax = plt.subplot(1, no_mscape, loop+1)


            if map_type == "perc":
                plot = ax.pcolormesh(samples, taxa, heatmap, vmin=0, vmax=100, cmap="magma_r")
            elif map_type == "hcid" :
                plot = ax.pcolormesh(samples, taxa, heatmap, vmin=min(df_min), vmax=max(df_max), cmap="Reds")
            else: #total count and top count
                plot = ax.pcolormesh(samples, taxa, heatmap, vmin=min(df_min), vmax=max(df_max), cmap="magma_r")
    

            #for hcid heatmap, add sample_id to x-axis, and total_counts to y-axis if there are contaminants in its column
            if map_type == "hcid":
                df = df.replace("NaN", 0)
                df = df.apply(convert_to_numeric)
                #only keep x-labels with counts
                xlabels = []
                for sample in df:
                    sample_count = df[sample].sum()
                    if int(sample_count) > 0:
                        xlabels.append(sample)
                    else:
                        xlabels.append("")
                ax.set_xticks(range(len(samples)),labels=xlabels, rotation=90)

                df["total"] = df.sum(axis=1)
                total = list(df["total"])
                total = [int(i) for i in total]
                hcid_list.append(total)

            else:
                plt.xticks([])

            # Minor ticks
            #ax.grid(color='w', linewidth=1.5)
            ax.set_xticks(np.arange(-.5, len(samples), 1), minor=True)
            ax.set_yticks(np.arange(-.5, len(taxa), 1), minor=True)
            # Gridlines based on minor ticks
            ax.grid(which='minor', color='w', linestyle='-', linewidth=1.5)

            if loop == 0:
                plt.yticks(taxa,rotation=0,fontsize='10')
                plt.xlabel(mscape_names[loop], fontweight='bold', horizontalalignment='center', rotation=90)

                if map_type != "hcid":
                    # Second X-axis
                    domain_dict = pd.DataFrame(domain_list, columns=["x"]).groupby('x').size().to_dict()

                    #Define the distances between each secondary y-axis tick
                    ytick_limits = [-0.5]
                    ytick_names = []
                    current_count = 0
                    ytick_loop = 0
                    for domain in domain_dict:
                        ytick_loop += 1
                        count = domain_dict[domain]

                        if map_type == "perc" or map_type == "top count":
                            ratio = (count * 0.97) + current_count
                        else:
                            if ytick_loop == len(domain_dict):
                                ratio = (count) + current_count - 0.5
                            else:
                                ratio = (count) + current_count

                        ytick_limits.append(ratio)
                        ytick_names.append(domain)
                        current_count = ratio


                    ytick_distance = []
                    last_point = 0
                    for limit in ytick_limits:
                        if limit == -0.5:
                            last_point = -0.5
                        else:
                            current_point = limit - last_point
                            halfway = current_point/2
                            label = last_point + halfway
                            ytick_distance.append(label)
                            last_point = limit
                
                    # label the classes:
                    sec = ax.secondary_yaxis(location=-3)
                    sec.set_yticks(ytick_distance, labels=ytick_names)
                    sec.tick_params('y', length=0, colors="k", grid_color = "k")

                    # lines between the classes:
                    sec2 = ax.secondary_yaxis(location=-3)
                    sec2.set_yticks(ytick_limits, labels=[], minor=False)
                    sec2.tick_params('y', which='major', length=10, width=1, direction="in", colors="k", grid_color = "k")

                    sec.spines["left"].set_visible(True)
                    sec.spines["left"].set_color("k")
                    sec.spines["left"].set_linewidth(1.5)
                    

                if map_type == "perc" or map_type == "top count":
                    ax.set_ylim(-0.5, 19.4)
                #else:
                    #ax.set_ylim(-0.5, len(genus)*0.97)

            elif loop < (no_mscape-1):
                plt.yticks([])      
                plt.xlabel(mscape_names[loop], fontweight='bold', horizontalalignment='center', rotation=90)
            elif loop == (no_mscape-1):
                total_hcid = [sum(x) for x in zip(*hcid_list)]

                if sum(total_hcid) > 0:
                    sec_y = ax.secondary_yaxis(location=1)
                    sec_y.set_yticks(range(len(taxa)), labels=total_hcid)
                    sec_y.tick_params('y', length=8)
                    hcid_display = "display: block"


                plt.yticks([])
                plt.xlabel(mscape_names[loop], fontweight='bold', horizontalalignment='center', rotation=90)
                cax = fig.add_axes([ax.get_position().x1+0.05,ax.get_position().y0,0.02,ax.get_position().height])
                cbar = plt.colorbar(plot, cax=cax)
                cbar.outline.set_color('black')
            loop += 1

        plt.subplots_adjust(wspace=0.02, hspace=0)
     

        # Save the figure as a Base64 string
        buf = io.BytesIO()
        plt.savefig(buf, format='png', bbox_inches='tight')
        buf.seek(0)
        heatmap = base64.b64encode(buf.read()).decode('utf-8')
        buf.close()
        plt.close(fig)

        return heatmap

    heatmap_list.append(plot_heatmaps(top_mscape_perc, "perc"))
    heatmap_list.append(plot_heatmaps(top_mscape_count, "top count"))
    heatmap_list.append(plot_heatmaps(total_mscape_count, "total count"))
    heatmap_list.append(plot_heatmaps(hcid_count, "hcid"))

    total_map_height = (total_mscape_count[0].shape[0])*22.5
    genera_count = total_mscape_count[0].shape[0]

    hcid_height = (hcid_count[0].shape[0])*22.5

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

                                absolute_class=stacked_class[0], relative_class=stacked_class[1],
                                top_perc_heatmap=heatmap_list[0], top_count_heatmap=heatmap_list[1], total_count_heatmap=heatmap_list[2], hcid_heatmap = heatmap_list[3], 
                                total_map_height=total_map_height, genera_count=genera_count, hcid_height=hcid_height, hcid_display=hcid_display)

    #os.makedirs(f'{output_path}/', exist_ok=True)
    print(output_path)
    # Save the rendered HTML to a file
    with open(f"{output_path}/summary_report.html", "w") as f:
        f.write(html_content)

    print(f"HTML file generated: {output_path}")

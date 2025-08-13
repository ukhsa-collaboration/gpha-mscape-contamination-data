#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import io
import base64
from mako.template import Template
import numpy as np
import seaborn as sns
import os
import argparse
import json
import scipy.stats as stats
from natsort import natsorted
from matplotlib.colors import LinearSegmentedColormap

from utils import define_heatmap_datasets, make_heatmap_df, make_date_list, convert_to_numeric, filter_multiples

def sort_by_date(count_df, date_list):
    #add and sort samples by date
    count_df.loc[len(count_df)] = date_list
    row_number = count_df.index.get_loc(count_df[count_df["Scientific_Name"] == "Date"].index[0])
    count_df.iloc[row_number, 2:-1] = pd.to_datetime(count_df.iloc[row_number, 2:-1], format='%Y-%m-%d')
    #Sort all samples by date
    count_df.iloc[row_number, 2:-1] = natsorted(count_df.iloc[row_number, 2:-1])

    #get timestamps line and turn to list
    timestamp_df = count_df.drop(columns="Counts_Overall")
    row_number = timestamp_df.index.get_loc(timestamp_df[timestamp_df["Scientific_Name"] == "Date"].index[0])
    date_df = timestamp_df.iloc[row_number]
    dates = list(date_df)
    #remove timestamps
    count_df = count_df[~count_df.Scientific_Name.str.contains("Date")]
    return count_df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse text files and a CSV file.")
    parser.add_argument('--reports', nargs='+', help="List of text files", required=True)
    parser.add_argument('--metadata', help="CSV file path", required=True)
    parser.add_argument('--site_key', help="JSON file specifying site name to number", required=True)
    parser.add_argument('--final_reports', help="Output directory", required=True)
    parser.add_argument('--template', help="HTMl template", required=True)
    parser.add_argument('--reference', help="excel file for dictionary of contaminants", required=True)
    #parser.add_argument('--hcids', nargs='+', help="List of CSV file for HCID contaminants", required=True)
    args = parser.parse_args()

    reports = args.reports
    metadata = pd.read_csv(args.metadata)
    template_dir = args.template
    #pip install openpyxl
    niche_table = pd.read_excel(args.reference)
    output_path = args.final_reports
    os.makedirs(output_path, exist_ok=True) #make directory path into real directory
    os.makedirs(f'{output_path}site_reports/', exist_ok=True)

    #site name to number
    site_key = {}
    with open(args.site_key, 'r') as f:
        site_key = json.load(f)

    grouped_metadata = metadata.groupby(['site', 'control_type_details'])

    datasets, mscape_datasets, samples, mscape_samples, sample_dates, mscape_dates = define_heatmap_datasets(grouped_metadata, site_key)

    exclude_cols = ["Scientific_Name"]

    present_dfs = []
    diff_loop = 0
    for needed_samples in mscape_samples:
        #for dataset in mscape_datasets:
        count_df, og_df = make_heatmap_df(needed_samples, reports, "count")

        sample_date_df = mscape_dates[diff_loop]
        date_list = make_date_list(sample_date_df)

        #name_line = ["domain", "scientific_name"] + sample_names[loop] + ["counts_overall"]
        #ßcount_df.loc[len(count_df)] = name_line
        #add and sort samples by date
        count_df.loc[len(count_df)] = date_list
        row_number = count_df.index.get_loc(count_df[count_df["Scientific_Name"] == "Date"].index[0])
        count_df.iloc[row_number, 2:-1] = pd.to_datetime(count_df.iloc[row_number, 2:-1], format='%Y-%m-%d')
        #Sort all samples by date
        count_df.iloc[row_number, 2:-1] = natsorted(count_df.iloc[row_number, 2:-1])

        #get timestamps line and turn to list
        timestamp_df = count_df.drop(columns="Counts_Overall")
        row_number = timestamp_df.index.get_loc(timestamp_df[timestamp_df["Scientific_Name"] == "Date"].index[0])
        date_df = timestamp_df.iloc[row_number]
        dates = list(date_df)
        #remove timestamps
        count_df = count_df[~count_df.Scientific_Name.str.contains("Date")]

        count_df = count_df.drop(columns=["Domain", "Counts_Overall"])

        include_cols = [col for col in count_df.columns if col not in exclude_cols]

        # Boolean indexing to filter rows showing only genus/family in the rank column
        filtered_df = count_df[count_df[include_cols].apply(lambda x: filter_multiples(x), axis=1)]

        present_dfs.append(filtered_df)
        diff_loop += 1



    #for each individual page
    site_loop = 0
    #for dataset in mscape_datasets:
    for needed_samples in mscape_samples:

        count_df, og_df = make_heatmap_df(needed_samples, reports, "count")
        print(count_df.columns)
        #count_df = count_df.drop(columns=["Taxon_ID"]) #drop taxon id column

        sample_date_df = mscape_dates[site_loop]
        date_list = make_date_list(sample_date_df)

        #name_line = ["domain", "scientific_name"] + sample_names[loop] + ["counts_overall"]
        #ßcount_df.loc[len(count_df)] = name_line
        
        count_df = sort_by_date(count_df, date_list)

        sorted_df = count_df.sort_values(by="Counts_Overall", ascending=False)

        def make_domain_df(sorted_df, domain):
            domain_df = sorted_df.loc[sorted_df["Domain"] == domain]

            top_df = domain_df.head(10) #top 10
            bottom_df = domain_df.tail(-10) #all taxa apart from top 10

            # Rearrange column 'Scientific Name' to the first position
            wordy_columns = ['Scientific_Name', 'Domain']

            # check if columns are not in the wordy_columns list or "counts_overall" (denoted by mscape_datasets)
            numeric_columns = [col for col in bottom_df.columns if col not in wordy_columns]
            top_df.loc["Other"] = bottom_df[numeric_columns].sum()
            top_df.at['Other', 'Domain'] = domain
            top_df.at['Other', 'Scientific_Name'] = "Other"

            top_df = top_df.drop(columns=["Domain", "Counts_Overall"])
            top_df.set_index("Scientific_Name", inplace=True)

            transposed_df = top_df.transpose()
            transposed_df.reset_index(inplace=True)

            # Turning it into percentages
            columns_to_sum = transposed_df.drop(columns=["index"])   
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
            relative_df = transposed_df.apply(make_relative)
            
            return transposed_df, relative_df

        absolute_dfs = []
        relative_dfs = []
        domains = ["Bacteria", "Fungi", "Archaea", "Viruses", "Protists"]
        for domain in domains:
            absolute_df, relative_df = make_domain_df(sorted_df, domain)
            absolute_dfs.append(absolute_df)
            relative_dfs.append(relative_df)

        no_sample = int(absolute_dfs[0].shape[0])

        def make_abundance_plots(domain_dfs):
            abundance_plots = []
            abundance_legends = []
            abundance_axes = []

            sns.set_style("whitegrid")
            domain_loop = 0
            for df in domain_dfs:
                # Generate the plot
                fig, ax = plt.subplots(figsize=(no_sample*0.8, 6))
                #hhttps://coolors.co/palette/001219-005f73-0a9396-94d2bd-e9d8a6-ee9b00-ca6702-bb3e03-ae2012-9b2226
                colour_list = ["#001219", "#005f73", "#0a9396", "#94d2bd", "#e9d8a6", "#ee9b00", "#ca6702", "#bb3e03", "#ae2012", "#9b2226"]

                to_be_removed = {'index'}
                genera = [item for item in list(df.columns) if item not in to_be_removed]
                colour_loop = 0
                colours = []
                for genus in genera:
                    if genus == "Other":
                        colours.append("#adb5bd")
                    else:
                        colours.append(colour_list[colour_loop])
                    colour_loop += 1
                
                df.plot(
                    x='index', kind='bar', stacked=True, ax=ax,
                    color=colours, width=0.98
                )
                plt.xticks([])
                plt.xlabel("Samples")
                plt.title(domains[domain_loop])
                legend = plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")

                # Save the legend as a Base64 string
                fig  = legend.figure
                fig.canvas.draw()
                bbox  = legend.get_window_extent()
                bbox = bbox.from_extents(*(bbox.extents + np.array([-5,-5,5,5])))
                bbox = bbox.transformed(fig.dpi_scale_trans.inverted())
                buf = io.BytesIO()
                plt.savefig(buf, format='png', bbox_inches=bbox)
                buf.seek(0)
                plot = base64.b64encode(buf.read()).decode('utf-8')
                abundance_legends.append(plot)
                buf.close()

                ax.get_legend().remove()
                # Save the figure as a Base64 string
                buf = io.BytesIO()
                plt.savefig(buf, format='png', bbox_inches='tight')
                buf.seek(0)
                plot = base64.b64encode(buf.read()).decode('utf-8')
                abundance_plots.append(plot)
                buf.close()

                # Save the y-axis as a Base64 string
                buf = io.BytesIO()
                plt.savefig(buf, format='png', bbox_inches=mtransforms.Bbox([[0, 0], [0.02, 1]]).transformed(
                        fig.transFigure - fig.dpi_scale_trans
                    ),
                )
                buf.seek(0)
                plot = base64.b64encode(buf.read()).decode('utf-8')
                abundance_axes.append(plot)
                buf.close()
                plt.close(fig)
                domain_loop += 1

            return abundance_plots, abundance_legends, abundance_axes

        abs_plots, abs_legends, abs_axes = make_abundance_plots(absolute_dfs)
        rel_plots, rel_legends, rel_axes = make_abundance_plots(relative_dfs)

        sns.set_style("dark")

        #organise mapped_required details for heatmap purposes


        #Label heatmap by pathogenicity
        sns.set_style("white")
        heatmaps = []

        #transpose niche_table data and subset for pathogenics
        patho_cols = ["Taxa", "Pathogenicity"]
        pathogenic_table = niche_table[patho_cols]
        pathogenic_table.fillna("NaN", inplace=True)
        pathogenic_table.set_index("Taxa", inplace=True)
        pathogenic_df = pathogenic_table.transpose()

        non_genus = ['root', 'unclassified']
        new_df = count_df[~count_df.Scientific_Name.isin(non_genus)]

        exclude_cols = ["Domain", "Niche", "Scientific_Name", "Counts_Overall"]
        include_cols = [col for col in new_df.columns if col not in exclude_cols]


        # Boolean indexing to filter rows showing only genus/family in the rank column
        new_df = new_df[new_df[include_cols].apply(lambda x: filter_multiples(x), axis=1)]

        pathogenicity = []
        for index, row in new_df.iterrows():
            genus_name = row["Scientific_Name"]
            
            if genus_name in pathogenic_df.columns:
                pathogen_str = list(pathogenic_df[genus_name])
                pathogenicity.extend(pathogen_str)
            else:
                pathogenicity.append("NaN")

        new_df["Pathogenicity"] = pathogenicity

        sorted_df = new_df.sort_values(by="Counts_Overall", ascending=False)
        sorted_df = sorted_df[(sorted_df["Counts_Overall"] > 0)]
        filtered_df = sorted_df.drop(columns=["Counts_Overall", "Domain"])

        sortfilt_df = filtered_df.sort_values(by="Pathogenicity", ascending=False)

        clean_df = sortfilt_df[~sortfilt_df.Pathogenicity.str.contains("NaN")]

        # get zeptometrix contamination from og_df (with all ranks)
        drop_df = og_df.drop(columns=["Taxon_ID", "Rank", "Domain"])
        og_date_list = date_list[1:-1]

        #add and sort samples by date
        drop_df.loc[len(drop_df)] = og_date_list
        row_number = drop_df.index.get_loc(drop_df[drop_df["Scientific_Name"] == "Date"].index[0])
        drop_df.iloc[row_number, 1:] = pd.to_datetime(drop_df.iloc[row_number, 1:], format='%Y-%m-%d')
        #Sort all samples by date
        drop_df.iloc[row_number, 1:] = natsorted(drop_df.iloc[row_number, 1:])

        #get timestamps line and turn to list
        row_number = drop_df.index.get_loc(drop_df[drop_df["Scientific_Name"] == "Date"].index[0])
        timestamps = drop_df.iloc[row_number]
        #remove timestamps
        drop_df = drop_df[~drop_df.Scientific_Name.str.contains("Date")]

        #get zeptometrix control df
        zepto = ["Bordetella", "Chlamydia", "Mycoplasmoides", "Adenoviridae", "Orthomyxoviridae", "Pneumoviridae", "Paramyxoviridae", "Picornaviridae", "Pneumoviridae", "Coronaviridae"]
        zepto_df = drop_df[drop_df.Scientific_Name.isin(zepto)]

        #sort zeptometrix df by total count
        num_zepto = zepto_df.drop(columns="Scientific_Name")
        zepto_df["total_count"] = num_zepto.sum(axis=1)
        zepto_df = zepto_df.sort_values(by="total_count", ascending=False)
        zepto_df = zepto_df.drop(columns="total_count")

        if len(zepto_df.index) == 0:
            empty_row = ["None"]
            loop = 1
            while loop < len(zepto_df.columns):
                empty_row.append(0)
                loop += 1
            zepto_df = pd.concat([pd.DataFrame([empty_row], columns=zepto_df.columns), zepto_df], ignore_index=True)

        #Get number of taxa 
        all_taxa = sortfilt_df.shape[0]
        nonan_taxa = clean_df.shape[0] #get the number of taxa recognised through lit review

        #remove zeptometrix taxa from clean_df
        clean_df = clean_df[~clean_df.Scientific_Name.isin(zepto)]

        df_list = []
        df_names = []

        grouped_df = clean_df.groupby(['Pathogenicity'])

        set_order = [["pathogenic", "potentially_pathogenic"], ["opportunistic", "potentially_opportunistic"], ["generally_commensal", "commensal"]]

        for set_type in set_order:
            current_df_list = []
            current_df_names = []
            for each_set in set_type:
                for sets in grouped_df:
                    df = sets[1]
                    name = list(sets[0])

                    if name[0] == each_set:
                        current_df_names.append(name[0])
                        df = df.drop(columns=["Pathogenicity"])
                        current_df_list.append(df)
            df_names.append(current_df_names)
            df_list.append(current_df_list)

        df_names = [["zeptometrix"]] + df_names
        df_list = [[zepto_df]] + df_list

        #make separate heatmaps for each subset of grouped_df
        annotation = []
        taxa_sums = []
        width = 0
        loop = 0
        for current_df_list in df_list:

            total_taxa = 0
            for df in current_df_list:
                taxa_no = df.shape[0]
                total_taxa = total_taxa + taxa_no
            taxa_sums.append(total_taxa)

            height_ratios = []
            for df in current_df_list:
                taxa_no = df.shape[0]
                height = taxa_no/total_taxa
                height_ratios.append(height)
            
            df_max = []
            df_min = []
            for df in current_df_list:
                print(df)
                df = df.drop(columns=["Scientific_Name"])
                df_max.append(max(df.max()))
                df_min.append(min(df.min()))

            #define figure dimensions
            height = total_taxa*0.3
            width = current_df_list[0].shape[1] * 0.32
            
            no_subplot = len(current_df_list)
            fig, ax = plt.subplots(figsize=(width,height), nrows=no_subplot, gridspec_kw={'height_ratios': height_ratios})

            patho_list = df_names[loop]

            name_loop = 0
            for map_df in current_df_list:

                map_df.set_index("Scientific_Name", inplace=True)
                map_df.replace(0, np.nan, inplace=True)
                
                map_df["total_counts"] = map_df.sum(axis=1)
                map_df = map_df.sort_values(by="total_counts", ascending=False)
                map_df = map_df.drop(columns="total_counts")

                genus = map_df.index.tolist()
                samples = map_df.columns.tolist()
                matrix = map_df.to_numpy()

                heatmap = np.reshape(matrix, (len(genus), len(samples)))
                heatmap = np.array(heatmap, dtype=np.float64)  # Ensure it's numeric
                
                ax = plt.subplot(no_subplot, 1, name_loop+1)

                current_patho = patho_list[name_loop]

                if current_patho == "zeptometrix":
                    custom_colors = ['#000000', '#5A1073']
                    shorthand = "Z"
                elif current_patho == "pathogenic":
                    # Custom color list for the colormap (e.g., shades of green, yellow, red)
                    custom_colors = ['#000000', '#ff0000']
                    shorthand = "P"
                elif current_patho == "potentially_pathogenic":
                    custom_colors = ['#000000', '#ff4d00']
                    shorthand = "PP"
                elif current_patho == "opportunistic":
                    custom_colors = ['#000000', '#ffbf00']
                    shorthand = "O"
                elif current_patho == "potentially_opportunistic":
                    custom_colors = ['#000000', '#ffde21']
                    shorthand = "PO"
                elif current_patho == "generally_commensal":
                    custom_colors = ['#000000', '#9dff00']
                    shorthand = "GC"
                elif current_patho == "commensal":
                    custom_colors = ['#000000', '#008000']
                    shorthand = "C"
                elif current_patho == "NaN":
                    custom_colors = ['#000000', '#0234d2']

                # Create a ListedColormap using the custom colors
                custom_cmap = LinearSegmentedColormap.from_list('pathogenicity', custom_colors)
                #plot = ax.pcolormesh(samples, genus, heatmap, cmap=custom_cmap, vmax=1000)
                im = ax.imshow(heatmap, vmin=min(df_min), vmax=max(df_max), cmap=custom_cmap)

                #plt.xticks([])
                #plt.yticks(genus,rotation=0,fontsize='10')
                #plt.ylabel(current_patho, labelpad=40)
                ax.set_xticks([])
                ax.set_yticks(range(len(genus)), labels=genus)

                #label the y-axis according to subplot dimensions
                if len(genus) >= 6:
                    ax.set_ylabel(current_patho, labelpad=40)
                else:
                    ax.set_ylabel(shorthand, labelpad=60)
                    annotation.append(f"{shorthand} = {current_patho}")

                cax = fig.add_axes([ax.get_position().x1+0.05,ax.get_position().y0,0.03,ax.get_position().height])
                cbar = plt.colorbar(im, cax=cax)
                
                #Get the max count in esch heatmap matrix
                num_map = np.nan_to_num(heatmap)
                max_count = num_map.max()

                new_labels = []
                label_loop = 0
                for element in cbar.ax.yaxis.get_ticklabels():
                    label_loop += 1
                    if label_loop == (len(cbar.ax.yaxis.get_ticklabels())):
                        max_count = int(max_count) + 1
                        new_labels.append(f"< {max_count}")
                    else:
                        label = element.get_text()
                        new_labels.append(element)
                
                cbar.ax.set_yticklabels(new_labels)

                name_loop += 1

            plt.subplots_adjust(wspace=0, hspace=0.01)
            buf = io.BytesIO()
            plt.savefig(buf, format='png', bbox_inches='tight')
            buf.seek(0)
            map = base64.b64encode(buf.read()).decode('utf-8')
            heatmaps.append(map)
            buf.close()
            plt.close(fig)

            loop += 1

        exclude_cols = ["Domain", "Niche", "Scientific_Name", "Counts_Overall"]
        include_cols = [col for col in count_df.columns if col not in exclude_cols]

        # Boolean indexing to filter rows showing only genus/family in the rank column
        filtered_df = count_df[count_df[include_cols].apply(lambda x: filter_multiples(x), axis=1)]

        nichemaps = []
        filtered_niche = filtered_df.drop(columns=["Domain"])
        colnames = list(filtered_niche["Scientific_Name"])

        filtered_table = niche_table.drop(columns=["Pathogenicity", "Reference"])
        filtered_table.set_index("Taxa", inplace=True)
        niche_df = filtered_table.transpose()

        all_niches = []
        for index, row in niche_df.iterrows():
            niche = list(niche_df.loc[index])
            new_niche = []
            for ele in niche:
                if isinstance(ele, str) and "," in ele:
                    ele = ele.split(",")
                    new_niche.extend(ele)
                elif isinstance(ele, str):
                    new_niche.append(ele)
            niche_list = list(set(new_niche))
            all_niches.append(niche_list)

        niche_colours = ['#77dd77','#ffc5d3', '#8fd3f4']
        
        niche_shorthand = {"soft_tissue": 'tissue', "resp_tract": 'resp',
                                "sterile_water": 'water', "taq_polymerase": 'taq', "lysing_enzymes": 'LE', "pcr_mix": 'PCR',
                                "extraction_kit": 'kit'}

        #create a dictionary for all niche categories, listing taxa they're found in as values
        niche_dict = {}
        count_dict = {}
        for taxon in niche_df:
            if taxon in colnames: #if taxon is in site, start looking for its niches from niche_df
                niche_list = list(niche_df[taxon])
                new_niche = []
                for niche in niche_list:
                    if isinstance(niche, str) and "," in niche:
                        niche = niche.split(",")
                        new_niche.extend(niche)
                    elif isinstance(niche, str):
                        new_niche.append(niche)
                #niche = [x for x in niche if str(x) != 'nan']
                for category in new_niche:
                    needed_line = filtered_niche[filtered_niche["Scientific_Name"] == taxon]
                    count = list(needed_line["Counts_Overall"])
                    if category in niche_dict:
                        niche_dict[category].append(taxon)
                        count_dict[category].extend(count)
                    else:
                        niche_dict[category] = [taxon]
                        count_dict[category] = count

        grouped_niches = []
        for niche_group in all_niches: #[lab], [human], [industry]
            site_by_niche = []
            for niche in niche_dict: #soil_water, primer, taq...
                current_niche = []
                if niche in niche_group: #e.g. soil_water in [soil_water, air...
                    for index, row in filtered_df.iterrows():
                        taxon_name = row["Scientific_Name"]
                        
                        if taxon_name in niche_dict[niche]: #if pseudomonas in niche_dict[soil_water...
                            current_niche.append(niche)
                        else:
                            current_niche.append("NaN")

                    filtered_df["Niche"] = current_niche
                    subset = filtered_df[filtered_df["Niche"] == niche]
                    subset = subset.sort_values(by="Counts_Overall", ascending=False)
                    subset = subset.drop(columns=["Domain"])
                    site_by_niche.append(subset)
            grouped_niches.append(site_by_niche)

        lab_taxa = []
        for every_niche in grouped_niches[0]: #find only niches in lab section
            taxa = list(every_niche["Scientific_Name"])
            lab_taxa.extend(taxa)

        lab_taxa = list(dict.fromkeys(lab_taxa))

        final_groups = [grouped_niches[0]]

        loop = 1
        while loop < len(grouped_niches):
            filtered_groups = []
            for group in grouped_niches[loop]:
                filtered_group = group.loc[~group["Scientific_Name"].astype(str).isin(lab_taxa)]
                if len(filtered_group.index) > 0:
                    filtered_groups.append(filtered_group)
            loop += 1
            final_groups.append(filtered_groups)

        #for every set of plots in the final grouping of niches
        niche_sums = []
        niche_annotations = []
        group_loop = 0
        for niche_group in final_groups:
            total_taxa = 0
            for df in niche_group:
                taxa_no = df.shape[0]
                total_taxa = total_taxa + taxa_no
            niche_sums.append(total_taxa)

            height_ratios = []
            for df in niche_group:
                taxa_no = df.shape[0]
                height = taxa_no/total_taxa
                height_ratios.append(height)
            
            df_max = []
            df_min = []
            for df in niche_group:
                df = df.drop(columns=["Scientific_Name", "Counts_Overall", "Niche"])
                df_max.append(max(df.max()))
                df_min.append(min(df.min()))

            height = total_taxa*0.3
            width = niche_group[0].shape[1] * 0.32
            
            no_subplot = len(niche_group)
            fig, axs = plt.subplots(figsize=(width,height), nrows=no_subplot, gridspec_kw={'height_ratios': height_ratios})

            custom_colors = ['#000000', niche_colours[group_loop]]

            df_loop = 0
            for map_df in niche_group:

                map_df.set_index("Scientific_Name", inplace=True)
                map_df.replace(0, np.nan, inplace=True)
                
                niche_name = list(map_df["Niche"])[0]
                map_df = map_df.drop(columns=["Niche", "Counts_Overall"])

                genus = map_df.index.tolist()
                samples = map_df.columns.tolist()
                matrix = map_df.to_numpy()
                
                heatmap = np.reshape(matrix, (len(genus), len(samples)))
                heatmap = np.array(heatmap, dtype=np.float64)  # Ensure it's numeric

                ax = plt.subplot(no_subplot, 1, df_loop+1)

                for niche in niche_colours:
                    if niche in niche_shorthand and len(genus) < 3: #if shorthand is needed
                        shorthand = niche_shorthand[niche]
                        niche_annotations.append(f'{shorthand} = {niche}')
       
                # Create a ListedColormap using the custom colors
                custom_cmap = LinearSegmentedColormap.from_list('niche', custom_colors)
                #plot = ax.pcolormesh(samples, genus, heatmap, cmap=custom_cmap)
                im = ax.imshow(heatmap, vmin=min(df_min), vmax=max(df_max), cmap=custom_cmap)

                #plt.xticks([])
                #plt.yticks(genus,rotation=0,fontsize='10')
                #plt.ylabel(current_patho, labelpad=40)
                ax.set_xticks([])
                ax.set_yticks(range(len(genus)), labels=genus)
                if niche_name in niche_shorthand and len(genus) < 3:
                    ax.set_ylabel(niche_shorthand[niche_name], labelpad = 40, rotation=0)
                else:
                    ax.set_ylabel(niche_name, labelpad = 40, rotation=0)
            
                if df_loop == (len(niche_group)-1):
                    #cax = fig.add_axes([ax.get_position().x1+0.1,ax.get_position().y0,0.02,ax.get_position().height])
                    cax = fig.add_axes([ax.get_position().x1+0.1, ax.get_position().y0+0.01, 0.02, 0.75])
                    cbar = plt.colorbar(im, cax=cax)


                df_loop += 1

            fig.align_ylabels(axs)
                              
            plt.subplots_adjust(wspace=0, hspace=0.02)
            buf = io.BytesIO()
            plt.savefig(buf, format='png', bbox_inches='tight')
            buf.seek(0)
            map = base64.b64encode(buf.read()).decode('utf-8')
            nichemaps.append(map)
            buf.close()
            plt.close(fig)
            group_loop += 1


        annotation = ', '.join(annotation)
        niche_annotations = ', '.join(niche_annotations)

        # Depict site-specific signals
        # Merge the DataFrames on a specific column
        all_df = pd.concat(present_dfs, axis = 0, join= "outer")  # Change join to 'outer' for outer join
        #merge on scientific name - this means that no scientific name is repeated if it's present in different dataframes
        all_df = all_df.groupby('Scientific_Name', as_index=False).first()
        all_df.fillna(0, inplace=True)

        new_dfs = []
        for df in present_dfs:
            new_df = all_df[df.columns]
            new_dfs.append(new_df)

        current_df = new_dfs[site_loop]
        other_dfs = [df for df in new_dfs if df is not current_df]

        current_name = mscape_datasets[site_loop]
        other_names = [name for name in mscape_datasets if name is not current_name]

        for df in other_dfs:
            df.set_index("Scientific_Name", inplace=True)
        current_df.set_index("Scientific_Name", inplace=True)

        all_taxa = list(current_df.index)
        for df in other_dfs:
            all_taxa.extend(df.index)

        all_taxa = list(set(all_taxa))

        p_dict = {}
        for taxa in all_taxa:
            if taxa in current_df.index: #only look for taxa that are in current_df
                taxa_p = []
                current_row = list(current_df.loc[taxa])
                for df in other_dfs:
                    if taxa in df.index:
                        comp_row = list(df.loc[taxa])
                    # print(np.var(current_row), np.var(comp_row))
                    # conduct the Welch's t-test
                    output = stats.ttest_ind(np.array(current_row), np.array(comp_row), equal_var = False)
                    output_list = list(output)
                    taxa_p.append(output_list[1])
                p_dict[taxa] = taxa_p

        p_values = pd.DataFrame(p_dict)

        p_index = []
        for name in other_names:
            p_index.append(f'{current_name} - {name}')

        p_values.index = p_index

        p_df = p_values.transpose()
        p_df = p_df.apply(convert_to_numeric).fillna(0).astype(float)

        sig_df = p_df[(p_df <= 0.05).any(axis=1)]

        site_name = mscape_datasets[site_loop]
        
        os.makedirs(f'{output_path}dataframes/', exist_ok=True)
        sig_df.to_csv(f'{output_path}/dataframes/{site_name}_ttest.csv', index=True)

        # Render the template with the Base64 string
        template = Template(filename=template_dir)
        html_content = template.render(site_name = site_name,
                                    bac_abs=abs_plots[0], fungi_abs=abs_plots[1], archaea_abs=abs_plots[2], virus_abs=abs_plots[3], protist_abs=abs_plots[4], no_sample=no_sample,
                                    bac_abun_l=abs_legends[0], fungi_abun_l=abs_legends[1], archaea_abun_l=abs_legends[2], virus_abun_l=abs_legends[3], protist_abun_l=abs_legends[4],
                                    bac_rel=rel_plots[0], fungi_rel=rel_plots[1], archaea_rel=rel_plots[2], virus_rel=rel_plots[3], protist_rel=rel_plots[4],
                                    zepto_map = heatmaps[0], patho_map = heatmaps[1], oppor_map = heatmaps[2], comme_map = heatmaps[3],
                                    zepto_count = taxa_sums[0], patho_count = taxa_sums[1], oppor_count = taxa_sums[2], comme_count = taxa_sums[3],
                                    all_taxa = all_taxa, nonan_taxa = nonan_taxa, width=width, annotation=annotation,
                                    lab_map = nichemaps[0], human_map = nichemaps[1], industry_map = nichemaps[2],
                                    lab_count = niche_sums[0],human_count = niche_sums[1], industry_count = niche_sums[2], niche_annotations=niche_annotations)

        # Save the rendered HTML to a file
        with open(f"{output_path}site_reports/{site_name}_report.html", "w") as f:
            f.write(html_content)

        print(f"HTML file generated: {output_path}/site_reports/")

        site_loop += 1


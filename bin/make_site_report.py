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
from natsort import natsorted
import json
import urllib.request
from matplotlib.colors import LinearSegmentedColormap

from utils import define_heatmap_datasets, make_heatmap_df, make_date_list

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse text files and a CSV file.")
    parser.add_argument('--reports', nargs='+', help="List of text files", required=True)
    parser.add_argument('--metadata', help="CSV file path", required=True)
    parser.add_argument('--site_key', help="JSON file specifying site name to number", required=True)
    parser.add_argument('--final_reports', help="Output directory", required=True)
    parser.add_argument('--template', help="HTMl template", required=True)
    parser.add_argument('--reference', help="CSV file for dictionary of contaminants", required=True)
    args = parser.parse_args()

    reports = args.reports
    metadata = pd.read_csv(args.metadata)
    template_dir = args.template
    niche_table = pd.read_csv(args.reference)
    output_path = args.final_reports
    os.makedirs(output_path, exist_ok=True) #make directory path into real directory
    os.makedirs(f'{output_path}site_reports/', exist_ok=True)

    #site name to number
    site_key = {}
    with open(args.site_key, 'r') as f:
        site_key = json.load(f)

    grouped_metadata = metadata.groupby(['site', 'control_type_details'])

    datasets, mscape_datasets, samples, mscape_samples, sample_dates, mscape_dates = define_heatmap_datasets(grouped_metadata, site_key)

    site_loop = 0
    #for dataset in mscape_datasets:
    for needed_samples in mscape_samples:

        count_df, og_df = make_heatmap_df(needed_samples, reports, "count")

        sample_date_df = mscape_dates[site_loop]
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

        sorted_df = count_df.sort_values(by="Counts_Overall", ascending=False)

        def make_domain_df(sorted_df, domain):
            domain_df = sorted_df.loc[sorted_df["Domain"] == domain]

            top_df = domain_df.head(10)
            bottom_df = domain_df.tail(-10)

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

        sns.set_style("darkgrid")

        with urllib.request.urlopen("https://raw.githubusercontent.com/artic-network/scylla/refs/heads/main/resources/hcid.json") as url:
            hcids = json.load(url)

        hcid_list = []
        for entry in hcids:
            name = entry["name"]
            name = name.replace(' ','_')
            hcid_list.append(name)

        hcid_df = og_df.loc[og_df['Scientific_Name'].isin(hcid_list)==True]

        hcid_df = hcid_df.drop(columns=["Domain", "Rank"])
        for hcid in hcid_list:
            hcid_exist = og_df.loc[og_df['Scientific_Name'].str.contains(hcid)]
            if len(hcid_exist.index) == 0:
                empty_hcid = [hcid]
                loop = 1
                while loop < len(hcid_df.columns):
                    empty_hcid.append(0)
                    loop += 1
                hcid_df = pd.concat([pd.DataFrame([empty_hcid], columns=hcid_df.columns), hcid_df], ignore_index=True)

        hcid_df.set_index('Scientific_Name', inplace=True)
        hcid_df.replace(0, np.nan, inplace=True)

        genus = hcid_df.index.tolist()
        samples = hcid_df.columns.tolist()
        matrix = hcid_df.to_numpy()

        fig, ax = plt.subplots(figsize=(18, 7))
        sns.heatmap(hcid_df, linewidth=1, cmap="OrRd") #cbar_kws={'label': 'Reads', 'orientation': 'horizontal'}
        hcidmap = np.reshape(matrix, (len(genus), len(samples)))
        hcidmap = np.array(hcidmap, dtype=np.float64)
        #plot = ax.pcolormesh(samples, genus, hcidmap, cmap='OrRd')


        #plt.colorbar(plot, ax=ax)
        plt.xticks([])
        plt.xlabel("Samples")
        plt.ylabel("")

        buf = io.BytesIO()
        plt.savefig(buf, format='png', bbox_inches='tight', pad_inches=0.1)
        buf.seek(0)
        hcidmap = base64.b64encode(buf.read()).decode('utf-8')
        buf.close()
        plt.close(fig)

        hcid_height = hcid_df.shape[0] * 23

        #Label heatmap by pathogenicity
        sns.set_style("white")
        heatmaps = []
        
        #transpose and rearrange table for data wrangling
        transposed_table = niche_table.transpose()
        transposed_table.reset_index(inplace=True)
        transposed_table.columns = transposed_table.iloc[0]
        transposed_table = transposed_table.drop(transposed_table.index[0])

        pathogenic_df = transposed_table.iloc[:1]
        pathogenic_df = pathogenic_df.dropna(axis=1, how='all')

        non_genus = ['root', 'unclassified']
        new_df = count_df[~count_df.Scientific_Name.isin(non_genus)]

        pathogenicity = []
        for index, row in new_df.iterrows():
            genus_name = row["Scientific_Name"]
            
            if genus_name in pathogenic_df.columns:
                pathogen_str = list(pathogenic_df[genus_name])
                pathogenicity.extend(pathogen_str)
            else:
                pathogenicity.append("NaN")

        new_df["Pathogenicity"] = pathogenicity
        #wordy_cols = ["Scientific_Name", "Pathogenicity"]
        #num_cols = new_df.drop(columns=wordy_cols)
        #num_cols["total_counts"] = num_cols.sum(axis=1)
        #new_df["total_counts"] = num_cols["total_counts"]

        sorted_df = new_df.sort_values(by="Counts_Overall", ascending=False)
        sorted_df = sorted_df[(sorted_df["Counts_Overall"] > 0)]
        filtered_df = sorted_df.drop(columns=["Counts_Overall", "Domain"])

        sortfilt_df = filtered_df.sort_values(by="Pathogenicity", ascending=False)

        clean_df = sortfilt_df[~sortfilt_df.Pathogenicity.str.contains("NaN")]

        #Get number of taxa 
        all_taxa = sortfilt_df.shape[0]
        nonan_taxa = clean_df.shape[0]

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

                if current_patho == "pathogenic":
                    # Custom color list for the colormap (e.g., shades of blue, green, red, and purple)
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
                im = ax.imshow(heatmap, cmap=custom_cmap)

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

        annotation = ', '.join(annotation)
        site_name = mscape_datasets[site_loop]

        # Render the template with the Base64 string
        template = Template(filename=template_dir)
        html_content = template.render(site_name = site_name,
                                    bac_abs=abs_plots[0], fungi_abs=abs_plots[1], archaea_abs=abs_plots[2], virus_abs=abs_plots[3], protist_abs=abs_plots[4], no_sample=no_sample,
                                    bac_abun_l=abs_legends[0], fungi_abun_l=abs_legends[1], archaea_abun_l=abs_legends[2], virus_abun_l=abs_legends[3], protist_abun_l=abs_legends[4],
                                    bac_rel=rel_plots[0], fungi_rel=rel_plots[1], archaea_rel=rel_plots[2], virus_rel=rel_plots[3], protist_rel=rel_plots[4],
                                    hcidmap=hcidmap, hcid_height=hcid_height,
                                    patho_map = heatmaps[0], oppor_map = heatmaps[1], comme_map = heatmaps[2],
                                    patho_count = taxa_sums[0], oppor_count = taxa_sums[1], comme_count = taxa_sums[2],
                                    all_taxa = all_taxa, nonan_taxa = nonan_taxa, width=width, annotation=annotation)

        # Save the rendered HTML to a file
        with open(f"{output_path}site_reports/{site_name}_report.html", "w") as f:
            f.write(html_content)

        print(f"HTML file generated: {output_path}/site_reports/")

        site_loop += 1


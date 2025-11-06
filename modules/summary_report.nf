#!/usr/bin/env nextflow

//take input as --reports and--metadata


/*
 * A python script which makes text files for R plots
 */
process make_shannon_script {
    label 'process_low'
    label "process_long"
    container 'community.wave.seqera.io/library/pip_numpy_pandas:426ad974eac1c1db'

    input:
    path reports
    path metadata
    path site_key

    output:
    path "text_files"

    script:
    """
    make_r_files.py --reports ${reports.join(' ')} --metadata ${metadata} --site_key ${site_key} --output_dir text_files/
    """
}


/*
 * An R script which produces output for shannon's diversity graphs
 */

process get_shannon_plot {

    label 'process_low'
    container 'community.wave.seqera.io/library/r-argparse_r-crayon_r-dplyr_r-ggplot2_r-vegan:eb552a73894bf74c'
    input:
    path text_files
    output:
    path "plots"

    script:
    """
    make_r_plots.R ${text_files}/ plots/
    """
}

/*
 * A process to rename the hcid files
 */
 
process renameHcids {

    container 'alpine:latest'
    
    stageAs( filePattern: String, value: Path )

    input: 
    tuple val(newName), path(s3file)
    
    output:
    path("${newName}"), emit: renamedHcids
    
        
    script:
    """
    echo "Work dir contents before rename:"
    ls
    #echo "Renaming ${s3file} ${newName}..."
    #mv ${s3file} ${newName}
    """
    
}    
    

/*
 * A python script which produces output for making the html report
 */

process make_summary_report {

    label 'process_low'
    container 'community.wave.seqera.io/library/pip_mako_matplotlib_natsort_pruned:44e99f335376fa3b'
    publishDir "${params.outdir}/", mode: 'copy' // Publish final report to local directory specified in params.config

    input:
    path reports
    path metadata
    path site_key
    path hcids
    path template
    path plots

    output:
    path "summary_reports/"

    script:
    """
    make_sum_report.py \
      --reports ${reports.join(' ')} \
      --metadata ${metadata} \
      --site_key ${site_key} \
      --plots_dir ${plots}/ \
      --hcids ${hcids.join(' ')} \
      --final_report summary_reports/ \
      --template ${template}
    """
}


workflow evaluate_negative_controls {
    report_list = params.reports?.split('\n') as List
    Channel
        .fromPath(report_list)
        .flatten()
        .collect()
        .set { reports }
    //reports.view()


    metadata_file = file(params.metadata, type: "file", checkIfExists:true)
    Channel
        .fromPath(metadata_file)
        .set { metadata }

    site_key = file(params.site_key, type: "file", checkIfExists:true)

    hcid_sample_list = params.hcids?.split('\n') as List
    println hcid_sample_list
    
    Channel
       .fromPath(hcid_sample_list, checkIfExists: true)
       .view { println "DEBUG: channel emits -> $it" }
       .map { file ->
           def climbid = file.parent.name
           //def filename = file.name
           //def newName = filename.contains(parent) ? filename : "${parent}.${filename}"
           tuple(climbid, file)
       }
       .set { hcidsToRename }
 
       //.flatten()
       //.collect()
       
    hcidsToRename.view { println "Input tuple: $it" }
       
    renameHcids(hcidsToRename)
    
    renameHcids.out.renamedHcids
        .collect()
        .set { all_renamed_hcid_files }
       
    template = file("$baseDir/bin/summary_report_template.html")

    make_shannon_script(reports, metadata, site_key)
    get_shannon_plot(make_shannon_script.out)
    make_summary_report(reports, metadata, site_key, all_renamed_hcid_files, template, get_shannon_plot.out)
    println "Report will be generated in ${params.outdir}"
}

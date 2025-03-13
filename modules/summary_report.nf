#!/usr/bin/env nextflow
params.script_path = "${workflow.projectDir}/bin"
params.metadata = "${params.metadata ?: ''}"
//take input as --reports and--metadata

/*
 * A Python script which parses the output of the previous script
 */
process make_shannon_script {

    container 'community.wave.seqera.io/library/pip_numpy_pandas:426ad974eac1c1db'

    input:
    path reports
    path metadata
    output:
    path "text_files/*"

    script:
    """
    python ${params.script_path}/make_r_files.py --reports ${reports.join(' ')} --metadata ${metadata} --output_dir text_files/
    """
}


/*
 * An R script which produces output for shannon's diversity graphs
 */

process get_shannon_plot {

    container 'community.wave.seqera.io/library/r-argparse_r-crayon_r-dplyr_r-ggplot2_r-vegan:eb552a73894bf74c'
    input:
    path "text_files/*"
    output:
    path "plots/*"

    script:
    """
    Rscript ${params.script_path}/make_r_plots.R text_files/ plots/
    """
}

/*
 * A python script which produces output for making the html report
 */

process make_report {

    container 'community.wave.seqera.io/library/pip_mako_matplotlib_natsort_pruned:44e99f335376fa3b'

    input:
    path reports
    path metadata
    path "plots/*"
    output:
    path "summary_report/*.html"

    publishDir "${params.output_dir}/", mode: 'copy' // Publish final report to local directory specified in params.config

    script:
    """
    python ${params.script_path}/make_sum_report.py --reports ${reports.join(' ')} --metadata ${metadata} --plots_dir plots/ --final_report summary_report/ --template ${params.script_path}/summary_report_template.html
    """
}


workflow evaluate_negative_controls {
    Channel
        .fromPath(params.reports)
        .collect()
        .set { reports }


    Channel
        .fromPath(params.metadata)
        .set { metadata }

    make_shannon_script(reports, metadata)
    get_shannon_plot(make_shannon_script.out)
    make_report(reports, metadata, get_shannon_plot.out)
    println "Report will be generated in ${params.output_dir}"
}

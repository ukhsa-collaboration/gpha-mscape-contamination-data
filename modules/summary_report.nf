#!/usr/bin/env nextflow

params.range = 100
params.script_path = "${workflow.projectDir}/../bin"

//take input as --input_dir
// Check if required parameters are provided
if (!params.input_dir) {
    exit 1, "Please provide --input_dir when running Nextflow."
}
/*
 * A Python script which parses the output of the previous script
 */
process make_shannon_script {
    input:
    path params.input_dir
    output:
    path "text_files/*"

    script:
    """
    python make_r_files.py ${params.input_dir} text_files/
    """
}

/*
 * An R script which produces output for shannon's diversity graphs
 */

process get_shannon_plot {
    input:
    path "text_files/*"
    output:
    path "plots/*"

    script:
    """
    Rscript make_r_plots.R text_files/ plots/
    """
}

/*
 * A python script which produces output for making the html report
 */

process make_report {
    input:
    path params.input_dir
    path "plots/*"
    output:
    path "report/*.html"

    publishDir "${System.getProperty('user.home')}/Downloads/", mode: 'copy', saveAs: {filename -> "negcontm_summary.html"} // Publish final report to local directory

    script:
    """
    python make_sum_report.py ${params.input_dir} plots/ report/ ${params.script_path}/summary_report_template.html
    """
}


workflow evaluate_negative_controls {
    make_shannon_script(params.input_dir)
    get_shannon_plot(make_shannon_script.out)
    make_report(params.input_dir, get_shannon_plot.out)
    println "Report will be generated in ~/Downloads/negcontm_summary.html"
}

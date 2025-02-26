#!/usr/bin/env nextflow

params.range = 100
//takes input as --input_dir
/*
 * A Python script which parses kraken files directory as input and produces taxa abundance tables
 */
process make_shannon_script {
    input:
    path params.input_dir
    output:
    path "text_files/*"

    script:
    """
    python ../bin/make_r_files.py ${params.input_dir} text_files/
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
    Rscript ../bin/make_r_plots.R text_files/ plots/
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

    publishDir ${file(params.input_dir).parent}, mode: 'copy', saveAs: {filename -> "summary_report.html"} // Publish final report to local directory

    script:
    """
    python ../bin/make_sum_report.py ${params.input_dir} plots/ report/
    """
}

workflow {
    make_shannon_script(params.input_dir)
    get_shannon_plot(make_shannon_script.out)
    make_report(params.input_dir, get_shannon_plot.out)
}

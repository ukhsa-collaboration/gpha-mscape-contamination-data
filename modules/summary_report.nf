#!/usr/bin/env nextflow

params.range = 100
//takes input as --input_dir
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
    python /Users/angelasun/Downloads/neg_control_reports/make_r_plots.py ${params.input_dir} text_files/
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
    Rscript /Users/angelasun/Downloads/neg_control_reports/shannon_plots.R text_files/ plots/
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

    publishDir '/Users/angelasun/Downloads/neg_control_reports/', mode: 'copy' // Publish final report to local directory

    script:
    """
    python /Users/angelasun/Downloads/neg_control_reports/make_summary_report.py ${params.input_dir} plots/ report/
    """
}

workflow {
    make_shannon_script(params.input_dir)
    get_shannon_plot(make_shannon_script.out)
    make_report(params.input_dir, get_shannon_plot.out)
}

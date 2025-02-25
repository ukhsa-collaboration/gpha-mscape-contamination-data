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
    python /https://raw.githubusercontent.com/aq-sun/contamination-data/main/bin/make_r_plots.py ${params.input_dir} text_files/
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
    Rscript /https://raw.githubusercontent.com/aq-sun/contamination-data/main/bin/shannon_plots.R text_files/ plots/
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

    publishDir '${params.input_dir}/*.html', mode: 'copy', saveAs: {filename -> "summary_report.html"} // Publish final report to local directory

    script:
    """
    python /https://raw.githubusercontent.com/aq-sun/contamination-data/main/bin/make_summary_report.py ${params.input_dir} plots/ report/
    """
}

workflow {
    make_shannon_script(params.input_dir)
    get_shannon_plot(make_shannon_script.out)
    make_report(params.input_dir, get_shannon_plot.out)
}

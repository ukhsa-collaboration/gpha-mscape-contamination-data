#!/usr/bin/env nextflow

params.script_path = "${workflow.projectDir}/../bin"

//take input as --kraken_reports
// Check if required parameters are provided
if (!params.kraken_reports) {
    exit 1, "Please provide --kraken_reports when running Nextflow."
}

Channel
    .fromPath(params.kraken_reports)
    .set { reports }

Channel
    .fromPath(params.metadata)
    .set { metadata }


/*
 * A Python script which parses the output of the previous script
 */
process make_shannon_script {
    input:
    path reports
    path metadata
    output:
    path "text_files/*"

    script:
    """
    python make_r_files.py --kraken_reports ${reports} --metadata ${metadata} --output_dir text_files/
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
    path reports
    path metadata
    path "plots/*"
    output:
    path "report/*.html"

    publishDir "${System.getProperty('user.home')}/Downloads/", mode: 'copy' // Publish final report to local directory

    script:
    """
    python make_sum_report.py --kraken_reports ${reports} --metadata ${metadata} --plots_dir plots/ --final_report report/ --template ${params.script_path}/summary_report_template.html
    """
}


workflow evaluate_negative_controls {
    make_shannon_script(reports, metadata)
    get_shannon_plot(make_shannon_script.out)
    make_report(reports, metadata, get_shannon_plot.out)
    println "Report will be generated in ~/Downloads/"
}

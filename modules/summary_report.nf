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
    python ${params.script_path}/make_r_files.py ${params.input_dir} text_files/
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
    Rscript ${params.script_path}/make_r_plots.R text_files/ plots/
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

    publishDir file(params.input_dir).getParent(), mode: 'copy', saveAs: {filename -> "summary_report.html"} // Publish final report to local directory

    script:
    """
    python ${params.script_path}/make_sum_report.py ${params.input_dir} plots/ report/ ${params.script_path}/summary_report_template.html
    """
}

workflow {
    make_shannon_script(params.input_dir)
    get_shannon_plot(make_shannon_script.out)
    make_report(params.input_dir, get_shannon_plot.out)

    // Automatically open the report after the workflow completes
    report_file = file("${params.input_dir}/../summary_report.html")
    if (report_file.exists()) {
        println "Report generated at: ${report_file}"
        "open ${report_file}".execute()  // For macOS
        "xdg-open ${report_file}".execute()  // For Linux
    }
}

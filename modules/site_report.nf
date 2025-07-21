#!/usr/bin/env nextflow

//take input as --reports and--metadata

/*
 * A python script which produces output for making the html report
 */

process rename {

    container 'community.wave.seqera.io/library/pip_mako_matplotlib_natsort_pruned:44e99f335376fa3b'

    input:
    tuple val(unique_id), path(hcid_counts)

    output:
    path "${unique_id}_hcid_counts.csv", emit: renamed_hcids

    script:
    """
    mv ${hcid_counts} "${unique_id}_hcid_counts.csv"
    """
}

process make_site_report {

    container 'community.wave.seqera.io/library/pip_mako_matplotlib_natsort_pruned:44e99f335376fa3b'

    input:
    path reports
    path metadata
    path hcids
    path site_key
    path template
    path reference

    output:
    path "site_reports/"

    publishDir "${params.outdir}/", mode: 'copy' // Publish final report to local directory specified in params.config

    script:
    """
    make_site_report.py \
      --reports ${reports.join(' ')} \
      --metadata ${metadata} \
      --hcids ${hcids.join(' ')} \
      --site_key ${site_key} \
      --final_report site_reports/ \
      --template ${template} \
      --reference ${reference}
    """
}


workflow evaluate_by_site {
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

    hcid_sample_list = params.hcids?.split('\n') as List
    Channel.from(hcid_sample_list)
    .map { full_path_str ->
        def file_path = file(full_path_str)
        def unique_id = file_path.getBaseName()
        def hcid_counts = file("${file_path}/hcid.counts.csv")
        return [unique_id, hcid_counts]
    }
    .set { hcid_sample_ch }

    hcid_sample_ch.view()
    rename(hcid_sample_ch)
    rename.out.renamed_hcids
       .flatten()
       .collect()
       .set { hcids }

    hcids.view()
    //hcids.view()
    //rename(hcids)
    

    site_key = file(params.site_key, type: "file", checkIfExists:true)

    template = file("$baseDir/bin/site_report_template.html")
    reference = file("$baseDir/bin/contaminant_dict.csv")

    make_site_report(reports, metadata, hcids, site_key, template, reference)
    println "Report will be generated in ${params.outdir}"
}

include { evaluate_by_site } from './modules/site_report'
include { evaluate_negative_controls } from './modules/summary_report'

// Check if required parameters are provided
if (!params.reports) {
    exit 1, "Please provide --reports when running Nextflow."
}

if (!params.metadata) {
    exit 1, "Please provide --metadata when running Nextflow."
}

if (!params.hcids) {
    exit 1, "Please provide --hcids when running Nextflow."
}

if (!params.site_key) {
    exit 1, "Please provide --site_key when running Nextflow."
}
workflow {
    evaluate_negative_controls()
    evaluate_by_site()
}
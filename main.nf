#!~/bin/nextflow

// params
params.dataset = "hartwig"
mode = params.dataset

params.index = "index.csv"

// prepare master channel
Channel
    .fromPath(params.index)
    .splitCsv(header:true)
    .map{ row-> tuple(row.sampleId,
                     file(row.vcfFile),
                     file(row.cnaFile),
                     row.purity) }
    .set { samples_ch }


process formatVCF{
    // label 'bcftools'
    publishDir 'tmp_files'
    
    input:
    set id, file(vcf), file(cna), val(purity) from samples_ch

    output:
    set id, file("${id}_pyclone_format.tsv"), val(purity) into  pyclone_format_ch

    script:
    if (mode == "hartwig")
        """
        bcftools query -i 'TYPE="snp"' \
                        -s ${id} \
                        -f '%CHROM:%POS:%REF:%ALT{0},%REF,%ALT{0},[%AD]\n' \
                        ${vcf} | parse_query_pcawg_sanger.awk > ${id}.tmpFile

        annotate_cna.R ${vcf} ${cna} ${id}.tmpFile hartwig ${id}

        rm ${id}.tmpFile
        """
    else if (mode == "pcawg_sanger")
        """
        bcftools query -s TUMOUR \
                       -f '%CHROM:%POS:%REF:%ALT{0} %REF %ALT{0} [ %FAZ %RAZ %FCZ %RCZ %FGZ %RGZ %FTZ %RTZ %PM %GT]\n' ${vcf} | parse_query_pcawg_sanger.awk > ${id}.tmpFile

        annotate_cna.R ${vcf} ${cna} ${id}.tmpFile sanger_pcawg TUMOUR

        rm ${id}.tmpFile
        """
     else
        error "Invalid dataset type: ${mode}"

}


pyclone_format_ch.into{pyclone_master_ch; stratification_master_ch }
process pyClone{
    publishDir "$params.outdir"
    label 'pyClone'

    input:
    set id, file(pyclone_file), val(purity) from pyclone_master_ch

    output:
    set id, file("${id}_pyClone_loci.tsv"), file("${id}_pyClone_clusters.tsv") into results_pyclone

    when:
    params.pyclone

    script:
    """
    # first step is to configure the analysis
    PyClone setup_analysis \
            --in_files ${pyclone_file} \
            --working_dir pyclone_wd \
            --tumour_contents ${purity} \
            --samples ${id}

    # second step is to run the analysis
    PyClone run_analysis --config_file pyclone_wd/config.yaml

    # get the summary table
    PyClone build_table --config_file pyclone_wd/config.yaml \
                        --out_file ${id}_pyClone_loci.tsv \
                        --table_type loci \
                        --max_clusters 4
    
    PyClone build_table --config_file pyclone_wd/config.yaml \
                        --out_file ${id}_pyClone_clusters.tsv \
                        --table_type cluster \
                        --max_clusters 4
    
    """

}

process computeStratification {
    publishDir "${params.outdir}/stratification"

    input:
    set id, file(input_file), val(purity) from stratification_master_ch

    output:
    set id, file("clonality_plots/${id}_clonality_ratio.pdf") into results_stratification_plots
    set id, file("clonality_results/${id}_mutations_strand_clonality.txt") into results_stratification

    when:
    params.stratification

    script:
    """
    clonality_single_sample.R -Pv -i $input_file -p $purity -s $id
    """
}






#!/usr/bin/env nextflow
/*

############  hyperclust ######################

hyperclust Pipeline.

by: David Mas-Ponte
aim: The pipeline aims to call mutation clusters in human tumors.
Documentation: https://github.com/davidmasp/hyperclust
*/

def helpMessage() {
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run davidmasp/hyperclust [OPT]
    Mandatory arguments:
      
      --index [file]          index file with vcf files and the corresponding samples.
      --genome [file]         masked fasta file to use as reference sequence.
      -profile [str]          configuration profile.
    Options:
      --dataset [str]         Different parsing options for common available datasets. (hartwig, pcawg_sanger and tcga_strelka)
      --stratification [bool] Wether or not to compute stratification for the samples
      --batchsize [int]       randommut specific option.
      --intraBS [int]         randommut specific option.
      --ws [int]              randommut specific option.
      --times [int]           randommut specific option.
      --clonal [bool]         clustmut specific option.
    Help:
      --help  [bool]         Show help
    """.stripIndent()
}
params.pairs = true
params.strand = true
params.clonal = true


params.help = false
// Show help message
if (params.help) {
    helpMessage()
    exit 0
}


/*
 * PARAMS
 */

params.dataset = "hartwig"
mode = params.dataset
params.index = "index.csv"
params.stratification = true

params.genome = "/g/strcombio/fsupek_data/users/dmas/data/ILUMINA_refseq/GRCh37/GRCh37_refseq_masked.fa"
params.assembly = "GRCh37"
genome = file(params.genome) // file here bc only 1!
params.batchsize = 10000
params.intraBS = 100
params.ws = 500000
params.times = 50

// todo ideally the combination of these channels should make all boosts
params.clonal = true

// prepare master channel
Channel
    .fromPath(params.index)
    .splitCsv(header:true)
    .map{ row-> tuple(row.sampleId,
                     file(row.vcfFile),
                     file(row.cnaFile),
                     row.purity) }
    .set{ samples_ch}


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
        # this is mixed vcf with indels
        # I should filter by PASS?
        # the M is max alleles and m is min alleles
        bcftools view -m2 -M2 -v snps -O z -o ${id}_snp.vcf.gz ${vcf}

        bcftools query -i 'TYPE="snp"' \
                        -s ${id} \
                        -f '%CHROM:%POS:%REF:%ALT{0},%REF,%ALT{0},[%AD]\n' \
                        ${id}_snp.vcf.gz | parse_query_hartwig.awk > ${id}.tmpFile


        annotate_cna.R ${id}_snp.vcf.gz ${cna} ${id}.tmpFile hartwig ${id}

        rm ${id}.tmpFile
        rm ${id}_snp.vcf.gz
        """
    else if (mode == "pcawg_sanger")
        """
        bcftools query -s TUMOUR \
                       -f '%CHROM:%POS:%REF:%ALT{0} %REF %ALT{0} [ %FAZ %RAZ %FCZ %RCZ %FGZ %RGZ %FTZ %RTZ %PM %GT]\n' ${vcf} | parse_query_pcawg_sanger.awk > ${id}.tmpFile

        annotate_cna.R ${vcf} ${cna} ${id}.tmpFile sanger_pcawg TUMOUR

        # this is not part of any channel, can be removed
        rm ${id}.tmpFile
        """
    else if (mode == "tcga_strelka")
        """
        bcftools view -m2 -M2 -v snps -O z -o ${id}_snp.vcf.gz ${vcf}

        bcftools query -s TUMOR \
                       -f '%CHROM:%POS:%REF:%ALT{0} %REF %ALT{0} %TQSS [ %DP %AU %CU %GU %TU]\n' \
                       ${id}_snp.vcf.gz | parse_query_TCGA_strelka.awk > ${id}.tmpFile

        annotate_cna.R ${id}_snp.vcf.gz ${cna} ${id}.tmpFile tcga_strelka TUMOR

        # this is not part of any channel, can be removed
        rm ${id}.tmpFile
        rm ${id}_snp.vcf.gz
        """
    else
        error "Invalid dataset type: ${mode}"

}

pyclone_format_ch.into{stratification_master_ch;  samples2_ch}

process computeStratification {
    publishDir "${params.outdir}/stratification"

    input:
    set id, file(input_file), val(purity) from stratification_master_ch

    output:
    set id, file("clonality_plots/${id}_clonality_ratio.pdf") into results_stratification_plots
    set id, file("clonality_results/${id}_mutations_strand_clonality.txt") into results_stratification_strand_clonality
    set id, file("clonality_results/${id}_mutations_pair_clonality.txt") into results_stratification_pair_clonality
    set id, file("clonality_results/${id}_mutations_clonality.txt") into results_stratification_clonality

    when:
    params.stratification

    script:
    """
    clonality_single_sample.R -Pv -i $input_file -p $purity -s $id
    """
}

// prepare input from the pyclone format
process formatRND{
    publishDir "TMP/"
    input: 
    set id, file(mutsFile), val(purity) from samples2_ch

    output:
    set id, file("${id}_rndFormat.tsv") into mutsfiles
    
    script:
    """
    cut -f 1 ${mutsFile} | tail -n +2 | \
        awk 'BEGIN {FS=":";OFS="\t"} {print(\$1,\$2,\$2,\$3,\$4,1,"$id")}' > \
        ${id}_rndFormat.tsv 

    """
}


// serialize the input genome
process serialize_genome {
    publishDir "TMP/"
    conda '/g/strcombio/fsupek_home/dmas/ENV/py36'
    afterScript 'set +u; conda deactivate'
    input:
    file "${params.assembly}.fa" from genome

    output:
    file "${params.assembly}.fa.p" into serial_genome
    file "available_chromosomes.txt" into available_chromosomes

    """
    python -m randommut -M serialize -g ${params.assembly}.fa -a ${params.assembly}

    egrep ">" ${params.assembly}.fa | sed 's/>//' > available_chromosomes.txt
    """
}


// Remove duplicated positions

process rmdup {
    publishDir 'TMP/'  

    tag "${big_file_dups}"

    input:
    set id, file(big_file_dups) from mutsfiles
    file available from available_chromosomes

    output:
    set id, file("${big_file_dups.baseName}_rmdup.tsv") into bigfiles_rmdup

    """
    egrep -f $available $big_file_dups | sort -k7,7 -k1,1 -k2n,2 -k5,5  | uniq > ${big_file_dups.baseName}_rmdup.tsv
    """
}

// Split input mut files
process split {
    publishDir 'TMP/'  

    tag "${big_file}"

    input:
    set id, file(big_file) from bigfiles_rmdup

    output:
    set id, file("${big_file}*") into split_files mode flatten

    """
    split -d -l ${params.batchsize} $big_file $big_file
    """
} // the second bigfile here is to work as prefix!

// check how the chanell looks like
// chanel of tuples with [id,splitfile]
// after we would need to groupTuple

// run split_files.subscribe{println it}


//Randomize splited file
process randomize {
    publishDir 'TMP/' 
    label 'py36'
    
    tag "${partial}"

    input:
    file genome from serial_genome
    set id, file(partial) from split_files

    output:
    set id, file("${partial}.randomized") into edited_files

    script:
    """
    python -m randommut -M randomize -g ${genome} -m ${partial} -a ${params.assembly} -o ${partial}.randomized -t ${params.times} -w ${params.ws} -b ${params.intraBS}
    """
}
 
//group back the channel
edited_files
  .groupTuple()
  .set { grouped_files }

// merge back output files 
process merge {
    publishDir "${params.outdir}/randommut"

    tag "${id}"

    input:
    set id, file(edited_files) from grouped_files

    output:
    set id, file("${id}_w${params.ws}.randomized.tsv") into randomized_files 

    script:
    """
    #!/bin/bash
    # obtained from https://stackoverflow.com/questions/24641948
    OutFileName=${id}_w${params.ws}.randomized.tsv
    i=0                                       # Reset a counter
    for filename in ${edited_files} ; do 
        if [ "\$filename"  != "\$OutFileName" ] ;      # Avoid recursion 
        then 
        if [[ \$i -eq 0 ]] ; then 
            head -1  \$filename >   \$OutFileName # Copy header 
        fi
        tail -n +2  \$filename >>  \$OutFileName # Append from the 2nd line
        i=\$(( \$i + 1 ))                        # Increase the counter
    fi
    done
    """
}

if( params.clonal ) {
    boosting_ch = results_stratification_strand_clonality
}
else {
    println "This option is still not implemented"
}

samples_boosting_ch = boosting_ch.join(randomized_files)

process clusterCallingBoost {
    tag "$sampleId - clustmut boost"
    publishDir "${params.outdir}/clustmut", mode: 'copy', overwrite: true
    // 123 is returned when no available mutations are found.
    validExitStatus 0,123
    input:
    set sampleId, file(boostingPath), file(samplePath) from samples_boosting_ch
    output:
    file "${sampleId}_strandClonality_distance_mutlist.txt" optional true into  cluster_calls_mlist_boost
    file "${sampleId}_strandClonality_distance_VRanges.rds" optional true into  cluster_calls_boost
    file "${sampleId}_strandClonality_plot.pdf" optional true 
    
    script:
    """
    clustmut distance -i . \
            --glob "*${samplePath}" \
            --recursive \
            -o "${sampleId}_strandClonality" \
            -N 1 \
            -b ${boostingPath} \
            -Vlv
    """
}

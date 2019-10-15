#!/usr/bin/awk -f

BEGIN {
    FS=","
    OFS="\t"
    print("mutation_id","ref_counts","var_counts","normal_cn")
}

{
    ref_reads=$4
    alt_reads=$5
    
    print($1,ref_reads,alt_reads,2)

}


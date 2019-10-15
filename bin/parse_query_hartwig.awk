#!/usr/bin/awk -f

BEGIN {
    OFS="\t"
    print("mutation_id","ref_counts","var_counts","normal_cn")
}

{   
    print($1,$4,$5,2)
}


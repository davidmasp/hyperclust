#!/usr/bin/awk -f

BEGIN {
    OFS="\t"
    print("mutation_id","ref_counts","var_counts","normal_cn")
}

{
    if ($2 == "A")
    {
        ref_reads=$4+$5
    } 
    else if ($2 == "C")
    {
        ref_reads=$6+$7
    }
    else if ($2 == "G")
    {
        ref_reads=$8+$9
    }
    else if ($2 == "T")
    {
        ref_reads=$10+$11
    }

    if ($3 == "A")
    {
        alt_reads=$4+$5
    } 
    else if ($3 == "C")
    {
        alt_reads=$6+$7
    }
    else if ($3 == "G")
    {
        alt_reads=$8+$9
    }
    else if ($3 == "T")
    {
        alt_reads=$10+$11
    }

    print($1,ref_reads,alt_reads,2)

}


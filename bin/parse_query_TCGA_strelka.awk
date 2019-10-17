#!/usr/bin/awk -f

# this is the query
# 1) '%CHROM:%POS:%REF:%ALT{0} 
# 2) %REF 
# 3) %ALT{0} 
# 4) %TQSS [ 
# 5) %DP 
# 6) %AU 
# 7) %CU 
# 8) %GU 
# 9) %TU]\n'

BEGIN {
    OFS="\t"
    print("mutation_id","ref_counts","var_counts","normal_cn")
}

{
    if ($2 == "A")
    {
        split($6,refA,",")
        ref_reads=refA[$4]
    } 
    else if ($2 == "C")
    {
        split($7,refC,",")
        ref_reads=refC[$4]
    }
    else if ($2 == "G")
    {
        split($8,refG,",")
        ref_reads=refG[$4]
    }
    else if ($2 == "T")
    {
        split($9,refT,",")
        ref_reads=refT[$4]
    }

    if ($3 == "A")
    {
        split($6,altA,",")
        alt_reads=altA[$4]
    } 
    else if ($3 == "C")
    {
        split($7,altC,",")
        alt_reads=altC[$4]
    }
    else if ($3 == "G")
    {
        split($8,altG,",")
        alt_reads=altG[$4]
    }
    else if ($3 == "T")
    {
        split($9,altT,",")
        alt_reads=altT[$4]
    }

    print($1,ref_reads,alt_reads,2)

}


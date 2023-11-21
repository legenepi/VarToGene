#!/usr/bin/awk -f

BEGIN {
    OFS = "\t"
}

NR == FNR {
    snp[$2] = $1
    next
}

#Chrom=$1
#Pos=$2
#Name=$3
#rsids=$4
#effectAllele=$5
#otherAllele=$6
#Beta=$7
#Pval=$8
#minus_log10_pval=$9
#SE=$10
#N=$11
#ImpMAF=$12

#Chrom	Pos	Name	rsids	effectAllele	otherAllele	Beta	Pval	minus_log10_pval	SE	N	ImpMAF
#chr1	23597	chr1:23597:G:A	NA	G	A	0.0766	0.872655	0.05916	0.477893	35287	0.00011

FNR == 1 {
    header = $0
    gene_id = gensub(/^.+\/[[:digit:]]+_[[:digit:]]+_(.+).b37.txt.gz/, "\\1", "g", INFILE)
}

FNR > 1 {
    id = gensub(/chr/, "", "g", $1)"_"$2"_"($5 < $6 ? $5 : $6)"_"($6 > $5 ? $6 : $5)
    if (id in snp && $8 < 0.05)
        result[id] = $0
}

END {
    if (length(result) > 0) {
        if (!LOOP || LOOP == 1)
            print "gene_id", "locus", "id", header
        for (i in result)
            print gene_id, snp[i], i, result[i]
    }
}

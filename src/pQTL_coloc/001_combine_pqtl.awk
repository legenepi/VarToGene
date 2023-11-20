#!/usr/bin/awk -f

#Rationale: combine all the UKB-pQTL with significant variant-pQTL association after look-up analysis in severe asthma GWAS
#credible set variants.

BEGIN {
    OFS = "\t"
}

NR == 1 {
    print "LOCUS", $0
}

FNR > 1 {
    sub("_ukb_protein_lookup.txt", "", FILENAME)
    print FILENAME, $0
}

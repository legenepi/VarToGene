#!/usr/bin/awk -f

# For each line of the first input file (the proxies) add a chr:pos key to the proxies array
# The array keys are unique so there won't be any duplicates even if you 
# add the same chr, pos twice.
#
# NR = line number over all input files
# FNR = line number in current file
# Hence when NR == FNR you must be in the first input file
# "next" means go to the next line in the file without processing any subsequent rules
NR == FNR {
    proxies[$2][$3] = 1
    next
}

# For the second input file (the pQTL results) output the header
# (processing will only get this far if the rule above didn't match)
FNR == 1 {
    print
}

# For subsequent lines in the the second file 
FNR > 1 {

    # Store a copy of the original line to output
    line = $0

    # Substitute ":" to split chr, pos into fields for matching
    gsub(":", "\t", $0)
}

# Print the original line if chr, pos is in our array of proxies
proxies[$1][$2] {
    print line
}

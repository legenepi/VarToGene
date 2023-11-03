#!/usr/bin/env Rscript

#Rationale: Use the locus2gene pipeline of OpenTarget to query 10 asthma studies published btw 2018-2023 and the fine-mapped
#sentinel variants found by our severe asthma GWAS. Output: L2G genes from these 10 studies.

library(plyr)
library(tidyverse)
library(purrr)
library(xlsx)
# Install relevant library for HTTP requests
library(httr)

# Set study_id and variant_id variables
study_id_list <- c("GCST010042","GCST90038616","FINNGEN_R6_J10_ASTHMA","GCST007798","SAIGE_495","GCST007995","GCST006911","GCST90018795","GCST009798","GCST90014325")
variant_id_list <- c("8_80380364_A_C")

#FUNCTION THAT RETURNS THE L2G GENES FOR EACH VARIANT-STUDY PAIR.

# Install relevant library for HTTP requests
library(httr)

# Set study_id and variant_id variables
study_id <- "GCST010042"
variant_id <- "8_80380364_A_C"

# Build query string
query_string = "
  query genePrioritisationUsingL2G($myVariantId: String!, $myStudyId: String! ){
    studyLocus2GeneTable(studyId: $myStudyId, variantId: $myVariantId){
      rows {
        gene {
          symbol
          id
        }
        yProbaModel
        yProbaDistance
        yProbaInteraction
        yProbaMolecularQTL
        yProbaPathogenicity
        hasColoc
        distanceToLocus
      }
    }
  }
"

# Set base URL of Genetics Portal GraphQL API endpoint
base_url <- "https://api.genetics.opentargets.org/graphql"

# Set variables object of arguments to be passed to endpoint
variables <- list("myStudyId" = study_id, "myVariantId" = variant_id)

# Construct POST request body object with query string and variables
post_body <- list(query = query_string, variables = variables)

# Perform POST request
r <- POST(url=base_url, body=post_body, encode='json')

# Print first entry of L2G data to RStudio console
head(content(r)$data$studyLocus2GeneTable$rows, 1)

lapply(content(r)$data$studyLocus2GeneTable$rows, function(x) x%>% select(gene))
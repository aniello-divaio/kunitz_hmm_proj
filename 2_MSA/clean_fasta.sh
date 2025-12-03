#!/bin/bash
awk '{if (substr($1,1,1)==">") {print "\n"toupper($1)} else {printf "%s",toupper($1)}}' kunitz_msa_22.fasta \
    | sed 's/PDB://g' \
    > kunitz_MSA_clean.ali

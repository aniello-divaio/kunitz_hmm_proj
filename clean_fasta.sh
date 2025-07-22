#!/bin/bash
awk '{if (substr($1,1,1)==">") {print "\n"toupper($1)} else {printf "%s",toupper($1)}}' efold_output.txt | sed 's/PDB://g' > kunitz_hmm_ready.fasta

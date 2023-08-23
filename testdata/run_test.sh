#!/usr/bin/bash
getorf -find 0 -table 1 -minsize 600 -sequence <(zcat metatranscriptomic_assembly.fna.gz) -outseq metatranscriptomic_assembly.faa
python ../iterativehmmsearch.py -q metatranscriptomic_assembly.getorf.faa -o test_rna_out -d rdrpscan_full.hmm --threads 24
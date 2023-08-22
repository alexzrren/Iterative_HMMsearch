rna_iterdir=/data/user/renzirui/Projects/WIV_ChinaBatVirome/VirusIdentification_1384/output_rna_v1
dna_iterdir=/data/user/renzirui/Projects/WIV_ChinaBatVirome/VirusIdentification_1384/output_dna_v1
wdir=/data/user/renzirui/Projects/WIV_ChinaBatVirome/VirusIdentification_1384/downstream_full
assembly_dir=/data/user/renzirui/Projects/WIV_ChinaBatVirome/VirusIdentification_1384/Assembly/

echo $iterdir

mkdir -p $wdir/sequences/putative_markerproteins
cat $rna_iterdir/Iteration_*/hmmsearch.domtblout | grep -v '#' | awk '{print $1}' | sort -u > $wdir/sequences/putative_markerproteins/putative_markerproteins.rna.list
cat $dna_iterdir/Iteration_*/hmmsearch.domtblout | grep -v '#' | awk '{print $1}' | sort -u > $wdir/sequences/putative_markerproteins/putative_markerproteins.dna.list
cat $wdir/sequences/putative_markerproteins/putative_markerproteins.rna.list $wdir/sequences/putative_markerproteins/putative_markerproteins.dna.list > $wdir/sequences/putative_markerproteins/putative_markerproteins.list
seqkit grep -f $wdir/sequences/putative_markerproteins/putative_markerproteins.list ${assembly_dir}*.faa >  $wdir/sequences/putative_markerproteins/putative_markerproteins.faa

mkdir -p $wdir/sequences/putative_viralcontigs
cat $wdir/sequences/putative_markerproteins/putative_markerproteins.list | grep -Po '\w+-k\d+_\d+' | sort -u > $wdir/sequences/putative_viralcontigs/putative_viralcontigs.list
cat $wdir/sequences/putative_markerproteins/putative_markerproteins.rna.list | grep -Po '\w+-k\d+_\d+' | sort -u > $wdir/sequences/putative_viralcontigs/rna.list
cat $wdir/sequences/putative_markerproteins/putative_markerproteins.dna.list | grep -Po '\w+-k\d+_\d+' | sort -u > $wdir/sequences/putative_viralcontigs/dna.list
seqkit grep -f $wdir/sequences/putative_viralcontigs/putative_viralcontigs.list ${assembly_dir}*.fna > $wdir/sequences/putative_viralcontigs/putative_viralcontigs.fna

echo "seqkit grep done!"

mkdir -p $wdir/hmmsearch_RdRpscan/
echo "hmmsearch -o /dev/null --noali --cpu 8 --domtblout hmmsearch_RdRpscan.domtblout /data/user/renzirui/Databases/RdRp-scan/full.hmm $wdir/sequences/putative_markerproteins/putative_markerproteins.faa" > $wdir/hmmsearch_RdRpscan/hmmsearch_rdrpscan.sh
qsub -clear -wd $wdir/hmmsearch_RdRpscan/ -q pathogendb.q -binding linear:8 -l vf=8g,num_proc=8 $wdir/hmmsearch_RdRpscan/hmmsearch_rdrpscan.sh

mkdir -p $wdir/blast/nr
echo "diamond blastx -q $wdir/sequences/putative_viralcontigs/putative_viralcontigs.fna -d /data/rsync/NR_2022-05.dmnd -o iterviralcontigs_nr.blastx --outfmt 6 qseqid qlen sseqid slen qstart qend sstart send evalue bitscore length pident mismatch gaps staxids sscinames sskingdoms skingdoms sphylums stitle qcovhsp scovhsp --evalue 1e-5 --block-size 32.0 --index-chunks 1 --threads 90" > $wdir/blast/nr/blastxNR.sh
qsub -clear -wd $wdir/blast/nr -q pathogendb.q -binding linear:45 -l vf=250g,num_proc=90 $wdir/blast/nr/blastxNR.sh

mkdir -p $wdir/blast/nt
echo "export BLASTDB=/jdfssz1/ST_HEALTH/P20Z10200N0206/kaixinyang/database/2022-05/NT/BLAST_INDEX
blastn -query $wdir/sequences/putative_viralcontigs/putative_viralcontigs.fna -out iterviralcontigs_nt.blastn -db /jdfssz1/ST_HEALTH/P20Z10200N0206/kaixinyang/database/2022-05/NT/BLAST_INDEX/nt -task megablast -evalue 1E-5 -max_target_seqs 50 -outfmt "'"6 qseqid qlen saccver slen qstart qend sstart send evalue bitscore length pident mismatch gaps staxid ssciname scomname sblastname sskingdom staxids sscinames scomnames sblastnames sskingdoms stitle qcovs qcovhsp qcovus"'" -num_threads 48" > $wdir/blast/nt/blastNT.sh
qsub -clear -wd $wdir/blast/nt -q pathogendb.q -binding linear:24 -l vf=200g,num_proc=48 $wdir/blast/nt/blastNT.sh

mkdir -p $wdir/blast/vhdb
echo "diamond blastx -o iterviralcontigs_vhdbcds.blastx -d /data/user/renzirui/Databases/Virushostdb/virushostdb.formatted.cds.faa.dmnd -q $wdir/sequences/putative_viralcontigs/putative_viralcontigs.fna --threads 12 --sensitive --max-target-seqs 25 --evalue 1E-5 --block-size 2.0 --index-chunks 1 --outfmt 6 qseqid qlen sseqid slen qstart qend sstart send evalue bitscore length pident mismatch gaps stitle qcovhsp scovhsp" > $wdir/blast/vhdb/blastx_vhdb.sh
qsub -clear -wd $wdir/blast/vhdb -q pathogendb.q -binding linear:12 -l vf=20g,num_proc=12 $wdir/blast/vhdb/blastx_vhdb.sh

echo "cd-hit-est -i $wdir/sequences/putative_viralcontigs/putative_viralcontigs.fna -o ANI90 -c 0.9 -aS 0.8 -g 1 -d 0 -T 8 -M 16384" > $wdir/sequences/putative_viralcontigs/cdhit_ani90.sh
#qsub -clear -wd $wdir/sequences/putative_viralcontigs -q pathogendb.q -binding linear:8 -l vf=16g,num_proc=8 $wdir/sequences/putative_viralcontigs/cdhit_ani90.sh

echo "cd-hit -i $wdir/sequences/putative_markerproteins/putative_markerproteins.faa -o AAI90 -c 0.9 -aS 0.8 -g 1 -d 0 -T 8 -M 16384" > $wdir/sequences/putative_markerproteins/cdhit_aai90.sh
#qsub -clear -wd $wdir/sequences/putative_markerproteins/ -q pathogendb.q -binding linear:8 -l vf=16g,num_proc=8 $wdir/sequences/putative_markerproteins/cdhit_aai90.sh

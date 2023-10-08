rna_iterdir=/ldfssz1/ST_HEALTH/P20Z10200N0206/P20Z10200N0206_pathogendb/zhaohailong/zhongda_SuddenDeath/baobi.virus/workdir/RNA_out
dna_iterdir=/ldfssz1/ST_HEALTH/P20Z10200N0206/P20Z10200N0206_pathogendb/zhaohailong/zhongda_SuddenDeath/baobi.virus/workdir/DNA_out
wdir=/jdfssz1/ST_HEALTH/P20Z10200N0206/renzirui/Analysis/SuddenDeath_test
assembly_dir=/jdfssz1/ST_HEALTH/P20Z10200N0206/renzirui/Analysis/SuddenDeath_test/

echo $iterdir

# mkdir -p $wdir/sequences/putative_markerproteins
# cat $rna_iterdir/Iteration_*/hmmsearch.domtblout | grep -v '#' | awk '{print $1}' | sort -u > $wdir/sequences/putative_markerproteins/putative_markerproteins.rna.list
# cat $dna_iterdir/Iteration_*/hmmsearch.domtblout | grep -v '#' | awk '{print $1}' | sort -u > $wdir/sequences/putative_markerproteins/putative_markerproteins.dna.list
# cat $wdir/sequences/putative_markerproteins/putative_markerproteins.rna.list $wdir/sequences/putative_markerproteins/putative_markerproteins.dna.list > $wdir/sequences/putative_markerproteins/putative_markerproteins.list
# seqkit grep -f $wdir/sequences/putative_markerproteins/putative_markerproteins.list ${assembly_dir}*.faa >  $wdir/sequences/putative_markerproteins/putative_markerproteins.faa

# mkdir -p $wdir/sequences/putative_viralcontigs
# cat $wdir/sequences/putative_markerproteins/putative_markerproteins.rna.list | perl -ne 'BEGIN{my %exist;} @a=split/_/; $out=join('_', @a[0..$#a-1]); if(!$exist{$out}){ print"$out\n"; $exist{$out}=1; }' > $wdir/sequences/putative_viralcontigs/rna.list
# cat $wdir/sequences/putative_markerproteins/putative_markerproteins.dna.list | perl -ne 'BEGIN{my %exist;} @a=split/_/; $out=join('_', @a[0..$#a-1]); if(!$exist{$out}){ print"$out\n"; $exist{$out}=1; }' > $wdir/sequences/putative_viralcontigs/dna.list
# cat $wdir/sequences/putative_viralcontigs/rna.list $wdir/sequences/putative_viralcontigs/dna.list > $wdir/sequences/putative_viralcontigs/putative_viralcontigs.list
# seqkit grep -f $wdir/sequences/putative_viralcontigs/putative_viralcontigs.list ${assembly_dir}*.fna > $wdir/sequences/putative_viralcontigs/putative_viralcontigs.fna

# echo "seqkit grep done!"

mkdir -p $wdir/hmmsearch_RdRpscan/
echo "hmmsearch -o /dev/null --noali --cpu 1 --domtblout hmmsearch_RdRpscan.domtblout /hwfssz5/ST_HEALTH/P20Z10200N0206/user/renzirui/Database/HMMs/RdRp/rdrpscan_full.hmm $wdir/sequences/putative_markerproteins/putative_markerproteins.faa" > $wdir/hmmsearch_RdRpscan/hmmsearch_rdrpscan.sh
qsub -clear -wd $wdir/hmmsearch_RdRpscan/ -P P20Z10200N0206 -q st.q -binding linear:1 -l vf=8g,num_proc=1 $wdir/hmmsearch_RdRpscan/hmmsearch_rdrpscan.sh

mkdir -p $wdir/blast/nr
echo "diamond blastx -q $wdir/sequences/putative_viralcontigs/putative_viralcontigs.fna -d /hwfssz5/ST_HEALTH/P20Z10200N0206/user/renzirui/Database/NCBI-BLAST/20230505/DIAMOND_NR/nr.dmnd -o iterviralcontigs_nr.blastx --outfmt 6 qseqid qlen sseqid slen qstart qend sstart send evalue bitscore length pident mismatch gaps staxids sscinames sskingdoms skingdoms sphylums stitle qcovhsp scovhsp --evalue 1e-5 --block-size 32.0 --index-chunks 1 --threads 8" > $wdir/blast/nr/blastxNR.sh
qsub -clear -wd $wdir/blast/nr -P P20Z10200N0206_super -q st_supermem.q -binding linear:4 -l vf=300g,num_proc=8 $wdir/blast/nr/blastxNR.sh

mkdir -p $wdir/blast/nt
echo "export BLASTDB=/hwfssz5/ST_HEALTH/P20Z10200N0206/user/renzirui/Database/NCBI-BLAST/20230505/NT
blastn -query $wdir/sequences/putative_viralcontigs/putative_viralcontigs.fna -out iterviralcontigs_nt.blastn -db /hwfssz5/ST_HEALTH/P20Z10200N0206/user/renzirui/Database/NCBI-BLAST/20230505/NT/nt -task megablast -evalue 1E-5 -max_target_seqs 50 -outfmt "'"6 qseqid qlen saccver slen qstart qend sstart send evalue bitscore length pident mismatch gaps staxid ssciname scomname sblastname sskingdom staxids sscinames scomnames sblastnames sskingdoms stitle qcovs qcovhsp qcovus"'" -num_threads 8" > $wdir/blast/nt/blastNT.sh
qsub -clear -wd $wdir/blast/nt -P P20Z10200N0206_super -q st_supermem.q -binding linear:4 -l vf=300g,num_proc=8 $wdir/blast/nt/blastNT.sh

mkdir -p $wdir/blast/vhdb
echo "diamond blastx -o iterviralcontigs_vhdbcds.blastx -d /hwfssz5/ST_HEALTH/P20Z10200N0206/user/renzirui/Database/VirusHostDB/release217/virushostdb.formatted.cds.dmnd -q $wdir/sequences/putative_viralcontigs/putative_viralcontigs.fna --threads 1 --sensitive --max-target-seqs 25 --evalue 1E-5 --block-size 2.0 --index-chunks 1 --outfmt 6 qseqid qlen sseqid slen qstart qend sstart send evalue bitscore length pident mismatch gaps stitle qcovhsp scovhsp" > $wdir/blast/vhdb/blastx_vhdb.sh
qsub -clear -wd $wdir/blast/vhdb -P P20Z10200N0206 -q st.q -binding linear:1 -l vf=16g,num_proc=1 $wdir/blast/vhdb/blastx_vhdb.sh

import pandas as pd
from Bio import SeqIO
import sys
import os
from datetime import datetime
from glob import glob
import time
import progressbar
sys.path.append("modules/")
import hmmtools as hmmtools
#import stockholmtools as stockholmtools
import seqtools as seqtools


def main():
    #Step1: hmm2seq
    hmmblocks = hmmtools.read_hmm(input='iter0_full.hmm')
    #print(hmmblocks)
    hmmtools.hmm2seq(hmmblocks)


def curtime():
        time = datetime.now().strftime('%T')+'   '
        return time


def runcommand(command):
    import subprocess
    p = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    stdout, stderr = p.communicate()
    code = p.returncode
    return(code, stdout, stderr)


def run_palmscan():
    (code, stdout, stderr) = runcommand("palmscan -search_pp %s -rdrp -ppout pp.fa -report pp.txt -fevout.pp.fev")
    return code


def parse_palmfev(fevfile):
    with open(fevfile, 'r') as fd:
        region_dict = dict()
        for line in fd.read().splitlines():
            linedict = dict(tuple(i.split('=')) for i in line.split())
            linedict['pp_region'] = linedict['pp_start']+':'+linedict['pp_end']
            region_dict[linedict['query']] = linedict
    return region_dict


def stockholm2blocks(stockholmfile):
    buffer = ''
    contentblocks = []
    with open(stockholmfile, 'r') as fd:
        for line in fd.readlines():
            if line.startswith('//'):
                buffer += line
                contentblocks.append(buffer)
                buffer = ''
            else:
                buffer += line
    sys.stderr.write(curtime()+'[INFO] Successfully read STOCKHOLM MSA of %d profiles\n' % (len(contentblocks)))
    return contentblocks   


def stockholm2fasta_perprofile(stockholmfile, outdir):
    count = 0
    for stockholmblock in stockholm2blocks(stockholmfile): 
        stockholmblockid = stockholmblock.splitlines()[1].split()[-1]
        stockholmblockfile = os.path.join(outdir, '.tmp_stockholm')
        with open(stockholmblockfile,'w') as fdout:
            fdout.write(stockholmblock)
        outfasta = os.path.join(outdir, stockholmblockid + '.aln')
        records = SeqIO.parse(stockholmblockfile, "stockholm")
        SeqIO.write(records, outfasta, "fasta")
        seqtools.msa2fasta(outfasta, os.path.join(outdir, stockholmblockid + '.fa'))
        count += 1
        
    os.remove(stockholmblockfile)
    return count

def hmmpress_cleanup(globpath):
    globresults = glob(globpath)
    if len(globresults):
        for old_hmmpress in glob(globpath):
            os.remove(old_hmmpress)
        sys.stderr.write(curtime()+'[WARN] Found old hmmpress results, removed!\n')


def cov_calc(range1, range2):
    #print(range1, range2)
    range1 = list(map(int, range1.split(':')))
    if range1[0] == range1[1]:
        return 0
    range2 = list(map(int, range2.split(':')))
    x = range(range1[0], range1[1])
    y = range(range2[0], range2[1])
    intersec_len = min(x[-1], y[-1])+1-max(x[0], y[0])
    cov = intersec_len/(range2[1]-range2[0])
    return cov


def rdrp_covfilter(domtblout, iter_time, coverage=0.8):
    dfhmm = pd.read_table(domtblout, comment='#', header=None, delim_whitespace=True)
    dfhmm.columns = 'target_name target_accession tlen query_name query_accession qlen seq_evalue seq_score seq_bias domain_# domain_of domain_c-evalue domain_i-evalue domain_score domain_bias hmm_start hmm_end aln_start aln_end env_start env_end acc desc'.split()
    dfhmm['hmm_coord'] = dfhmm['hmm_start'].astype(str) + ':' + dfhmm['hmm_end'].astype(str)
    if iter_time > 1:
        dfhmm['query_name'] = dfhmm['query_name'].apply(lambda x: x.split('fasta_')[0]+'fasta')
    df_filtered = pd.DataFrame()
    for query_name, dfquery in dfhmm.groupby('query_name'):
        try:
            query_range = dfquery[dfquery['target_name']==query_name].reset_index().loc[0,'hmm_coord']
        except:
            continue
        dfquery['cov'] = dfquery['hmm_coord'].apply(lambda x: cov_calc(x, query_range))
        filtered_dfquery = dfquery[dfquery['cov']>=coverage]
        #print(filtered_dfquery)
        df_filtered = pd.concat([df_filtered, filtered_dfquery])
    return df_filtered


def seqkit_grep(seqid_list, fastain, fastaout):
    with open('.tmpseqid.list', 'w') as fdout:
        fdout.writelines(list(map(lambda x: x+'\n',seqid_list)))
    (code, _, stderr) = runcommand("cat %s | sed 's/\/.*$//g' | seqkit grep -f .tmpseqid.list > %s" % (fastain, fastaout))
    if code:
        sys.stderr.write(curtime()+'[ERRO] seqkit exit with code %d\n' % code)
        sys.stderr.write(stderr.decode())
        sys.exit(1)


def iterative_hmmsearch(query, ncpu=4, wdir='/home/renzirui/Analysis/Iterative_hmmsearch'):
    sys.stderr.write('-----------------------Iterative_hmmsearch-------------------------\n')
    if not os.path.exists(os.path.join(wdir, 'builded_hmm')):
        os.mkdir(os.path.join(wdir, 'builded_hmm'))
    for i in range(1, 11):
        sys.stderr.write('##Iteration %d\n' % i)
        iter_wdir = os.path.join(wdir, 'Iteration_%d' % i)
        if not os.path.exists(iter_wdir):
            os.mkdir(iter_wdir)
        else:
            sys.stderr.write(curtime()+'[WARN] Path Iteration_%d/ already exists! Files will be overwritten!\n' % i)

        #Step00. hmmpress
        _t0 = time.time()
        hmmpress_cleanup("{wdir}/builded_hmm/iter{last_iter}_full.hmm.*".format(wdir=wdir, last_iter=str(i-1)))
        code = os.system('/home/renzirui/miniconda3/bin/hmmpress {wdir}/builded_hmm/iter{last_iter}_full.hmm > /dev/null 2>&1'.format(wdir=wdir, last_iter=str(i-1)))
        if code:
            sys.stderr.write(curtime()+'[ERRO] hmmpress exit with code %d\n'%code)
        _t1 = time.time()
        sys.stderr.write(curtime()+'[INFO] hmmpress done (Time Elapsed: %.3fs)\n' % (_t1-_t0))
        sys.stderr.write(curtime()+'[INFO] HMM search iteration %d running...\n' % i)

        #Step01. hmmsearch_run
        _t2 = time.time()
        os.system('/home/renzirui/miniconda3/bin/hmmsearch --cpu {ncpu} -E 1e-10 -o /dev/null -A {wdir}/Iteration_{iter}/hmmsearch.aln --tblout {wdir}/Iteration_{iter}/hmmsearch.tblout --domtblout {wdir}/Iteration_{iter}/hmmsearch.domtblout {wdir}/builded_hmm/iter{last_iter}_full.hmm {query}'.format(ncpu=ncpu, wdir=wdir, query=query, iter=i, last_iter=i-1))
        _t3 = time.time()
        (code, stdout, _) = runcommand("cat %s/hmmsearch.domtblout | grep -v '#' | awk '{print $1}' | sort -u | wc -l" % (iter_wdir))
        num_hitseq = stdout.decode().strip()
        sys.stderr.write(curtime()+'[INFO] hmmsearch done, %s sequences hit (Time Elapsed: %.3fs)\n' % (num_hitseq, _t3-_t2))

        #Step02. Confirm existence of RdRp Motif
        #Step02.1 STOCKHOLM to FASTA (per profile)
        sto2fas_path = os.path.join(iter_wdir, 'sto2fas')
        try:
            os.mkdir(sto2fas_path)
        except FileExistsError:
            sys.stderr.write(curtime()+'[WARN] path %s exists, content will be overwritten!\n' % sto2fas_path)
        count = stockholm2fasta_perprofile("{iter_wdir}/hmmsearch.aln".format(iter_wdir=iter_wdir), sto2fas_path)
        _t4 = time.time()
        sys.stderr.write(curtime()+"[INFO] Converted %d records from STOCKHOLM format to FASTA (Time Elapsed: %.3fs)\n" % (count, _t4-_t3))

        ##Step02.2 Confirm RdRp Motifs(>80% cov sequence will be kept)
        filtered_df = rdrp_covfilter(iter_wdir+'/hmmsearch.domtblout', i, coverage=0.8)
        #print(filtered_df)
        num_filtered = len(filtered_df['target_name'].drop_duplicates())
        _t5 = time.time()
        sys.stderr.write(curtime()+"[INFO] %d sequences covered >80%% RdRp motif region (Time Elapsed: %.3fs)\n" % (num_filtered, _t5-_t4))

        ##Step02.3 Fetch Confirmed Sequences
        num_profile, num_seq = 0, 0
        for profile, profile_df in filtered_df.groupby('query_name'):
            num_profile += 1
            seqlist = list(profile_df['target_name'].drop_duplicates())
            num_seq += len(seqlist)
            if i > 1:
                seqkit_grep(seqlist, os.path.join(iter_wdir, 'sto2fas', profile+'_iter%d.fa'%(i-1)), os.path.join(iter_wdir, 'sto2fas', profile+'.filt.fa'))
            else:
                seqkit_grep(seqlist, os.path.join(iter_wdir, 'sto2fas', profile+'.fa'), os.path.join(iter_wdir, 'sto2fas', profile+'.filt.fa'))
        _t6 = time.time()
        sys.stderr.write(curtime()+"[INFO] %d sequences from %d profiles after filtering written to FASTA (Time Elapsed: %.3fs)\n" % (num_seq, num_profile, _t6-_t5))

        ##Step03.1 CD-HIT Clustering Sequences
        cdhit_path = os.path.join(iter_wdir, 'cluster')
        try:
            os.mkdir(cdhit_path)
        except FileExistsError:
            sys.stderr.write(curtime()+'[WARN] path %s exists, content will be overwritten!\n' % cdhit_path)
        bar = progressbar.ProgressBar(widgets=[curtime(), '[INFO] cd-hit clustering profile ',progressbar.SimpleProgress(), ' (', progressbar.Percentage(), ')'])
        for file in bar(glob(iter_wdir + '/sto2fas/*.filt.fa')):
            fileprefix = '.'.join(file.split('/')[-1].split('.')[:-2])
            code = os.system('cd-hit -d 0 -c 0.8 -aS 0.5 -i {query} -o {result_prefix} -T 4 > {result_prefix}.cdhit.log 2>&1'.format(query=file, result_prefix=os.path.join(cdhit_path, fileprefix)))
            if code:
                sys.stderr.write(curtime()+'[ERRO] cd-hit exit with code %s\n' % code)
        _t7 = time.time()
        sys.stderr.write(curtime()+'[INFO] cd-hit done (Time Elapsed: %.3fs)\n' % (_t7-_t6))

        ##Step03.2 Fetch Clustered Sequences
        hmmbuild_msa_path = os.path.join(iter_wdir, 'hmmbuild_MSA')
        try:
            os.mkdir(hmmbuild_msa_path)
        except FileExistsError:
            sys.stderr.write(curtime()+'[WARN] path %s exists, content will be overwritten!\n' % hmmbuild_msa_path)
        _num_seq = 0
        for _count, file in enumerate(glob(cdhit_path + '/*.clstr')):
            dfclus = seqtools.read_cdhit(file, '0.8')
            profile = '.'.join(file.split('/')[-1].split('.')[:-1])
            keepseqid = list(dfclus.loc[dfclus['Cluster_Sign_0.8']=='*','Sequence_ID'])
            _num_seq += len(keepseqid)
            if i > 1:
                seqkit_grep(keepseqid, os.path.join(iter_wdir,'sto2fas', profile+'_iter%d.aln'%(i-1)), os.path.join(hmmbuild_msa_path, profile+'_iter%d.aln'%i))
            else:
                seqkit_grep(keepseqid, os.path.join(iter_wdir,'sto2fas', profile+'.aln'), os.path.join(hmmbuild_msa_path, profile+'_iter%d.aln'%i))

        _t8 = time.time()
        sys.stderr.write(curtime()+'[INFO] %d MSA seqs for %d profiles fetched (Time Elapsed: %.3fs)\n' % (_num_seq, _count+1, _t8-_t7))

        ##Step04. hmmbuild
        hmmbuild_path = os.path.join(iter_wdir, 'hmmbuild')
        try:
            os.mkdir(hmmbuild_path)
        except FileExistsError:
            sys.stderr.write(curtime()+'[WARN] path %s exists, content will be overwritten!\n' % hmmbuild_path)
        bar2 = progressbar.ProgressBar(widgets=[curtime(), '[INFO] hmmbuild running on profile ',progressbar.SimpleProgress(), ' (', progressbar.Percentage(), ')'])
        for file in bar2(glob(hmmbuild_msa_path+'/*')):
            fileprefix = file.split('/')[-1]
            code = os.system('hmmbuild {hmm_prefix}.hmm {msa_prefix} > {hmm_prefix}.log 2>&1'.format(hmm_prefix=os.path.join(hmmbuild_path, fileprefix), msa_prefix=os.path.join(hmmbuild_msa_path, fileprefix)))
            if code:
                sys.stderr.write(curtime()+'[ERRO] hmmbuild exit with code %s\n' % code)

        code = os.system('cat {hmmbuild_path}/*.hmm > {iterative_hmm}'.format(hmmbuild_path=hmmbuild_path, iterative_hmm=os.path.join(wdir, 'builded_hmm', 'iter%d_full.hmm'%i)))
        sys.stderr.write(curtime()+'[INFO] Full iterative hmm file written to %s\n' % (os.path.join(wdir, 'builded_hmm', 'iter%d_full.hmm'%i)))
        _t9 = time.time()
        sys.stderr.write(curtime()+'[INFO] hmmbuild done (Time Elapsed: %.3fs)\n' % (_t9-_t8))
        sys.stderr.write(curtime()+'[INFO] Iteration 1 finished (Time Elapsed: %.3fs)\n\n\n' % (_t9-_t1))
        

if __name__ == '__main__':
    #main()
    iterative_hmmsearch('/home/renzirui/Analysis/Iterative_hmmsearch/allORFs_clean_withPublic.faa', ncpu=12, wdir='/home/renzirui/Analysis/Iterative_hmmsearch_v2')
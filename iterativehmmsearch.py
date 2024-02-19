import sys
import os
import datetime
import time
import glob
import time
import pandas as pd
from Bio import SeqIO
from io import StringIO
import argparse
import subprocess
# import custom scripts
sys.path.append(os.path.join(os.path.split(sys.argv[0])[0], "modules/"))
import progressbar
import hmmtools
import seqtools


def main():
    parser = cmd_argparse()
    args = parser.parse_args()
    if not argcheck(args):
        parser.print_usage(file=sys.stderr)
        sys.exit(1)
    if not prerequisites():
        sys.exit(1)
    printargs(args)
    cleanquery = preprocess(args.query, args.output, args.hmmdb)
    iterative_hmmsearch(cleanquery, args.output, args.threads, args.iteration, args.evalue, args.motifcov)
    

def argcheck(args):
    if not os.path.exists(os.path.join(os.getcwd(), args.query)):
        sys.stderr.write('[ERRO] Cannot access query input %s: No such file or directory\n' % (args.query))
        return False
    else:
        if not os.access(args.query, os.R_OK):
            sys.stderr.write('[ERRO] Cannot access query input %s: Permission denied\n' % (args.query))
            return False
            
    if not os.path.exists(args.output):
        try:
            os.mkdir(args.output)
            sys.stderr.write('[INFO] Output directory created\n')
        except:
            if not os.path.isdir(args.output):
                sys.stderr.write('[ERRO] Cannot write to output directory %s: Is a File, not a directory\n' % (args.output))
                return False
            elif not os.access(args.output, os.W_OK):
                sys.stderr.write('[ERRO] Cannot write to output directory %s: Permission denied\n' % (args.output))
                return False
            else:
                sys.stderr.write('[ERRO] Cannot access output directory %s: No such file or directory\n' % (args.output))
                return False
        
    if not os.path.exists(os.path.join(os.getcwd(), args.hmmdb)):
        sys.stderr.write('[ERRO] Cannot access hmmdb %s: No such file or directory\n' % (args.hmmdb))
        return False
    else:
        if not os.access(args.hmmdb, os.R_OK):
            sys.stderr.write('[ERRO] Cannot access hmmdb %s: Permission denied\n' % (args.hmmdb))
            return False
    return True
    

def cmd_argparse():
    parser = argparse.ArgumentParser()
    parser.add_argument('-q', '--query', metavar='<FASTA>', help='input query orf sequences in FASTA format', required=True)
    parser.add_argument('-o', '--output', metavar='<OUTPUT_DIRECTORY>', help='output and intermediate file directory', required=True)
    parser.add_argument('-d', '--hmmdb', metavar='<HMM_FILE>', help='hmm profile of ', required=True)
    parser.add_argument('--threads', help='threads number assigned in iterative hmmsearch [4]', default=4, type=int, metavar='INT')
    parser.add_argument('--evalue', metavar='FLOAT', help='Evalue threshold in iterative hmmsearch [1e-10]', default='1e-10', type=str)
    parser.add_argument('--iteration', metavar='INT', help='Iteration times of hmmsearch [10]', default=10, type=int)
    parser.add_argument('--motifcov', metavar='FLOAT', help='Motif coverage cutoff of hmmsearch results [0.75]', default=0.75, type=float)
    parser.add_argument('--iterstart', metavar='INT', help='Manually assign start point of iteration for additional iteration or rerun at breakpoint [1]', default=1, type=int)
    return parser


def check_command(command):
    code, _ = subprocess.getstatusoutput(command)
    if code == 127:
        return False
    else:
        return True


def prerequisites():
    command_used = ['seqkit', 'hmmsearch', 'hmmbuild', 'hmmpress', 'cd-hit']
    for command in command_used:
        if not check_command(command):
            print('[ERRO] %s not found, please install it first or add it into $PATH manually' % command)
            return False
    return True


def printargs(args):
    print('IterativeHMM_Searcher v0.2a\n\nAuthor:Zirui Ren <renzirui@genomics.cn>\n\n###  alpha version  ###\n')
    print('='*35+'\n'+"       Search Configurations\n"+'-'*35)
    print("%+13s%s%s" % ('Query',' : ', args.query))
    print("%+13s%s%s" %("Output"," : ", args.output))
    print("%+13s%s%s" %("HMMdb"," : ", args.hmmdb))
    print("%+13s%s%s" %("Threads"," : ", args.threads))
    print("%+13s%s%s" %("IterationNum"," : ", args.iteration))
    print("%+13s%s%s" %("MotifCov"," : ", args.motifcov))
    print("%+13s%s%s\n" %("IterStart"," : ", args.iterstart)+ '='*35 + '\n')


def curtime():
        time = datetime.datetime.now().strftime('%T')+'   '
        return time


def runcommand(command):
    import subprocess
    p = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    stdout, stderr = p.communicate()
    code = p.returncode
    return(code, stdout, stderr)


def parse_palmfev(fevfile):
    with open(fevfile, 'r') as fd:
        region_dict = dict()
        for line in fd.read().splitlines():
            linedict = dict(tuple(i.split('=')) for i in line.split())
            linedict['pp_region'] = linedict['pp_start'] + ':' + linedict['pp_end']
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
        massive = False
        if stockholmblockid == '"<seq#>|"':
            stockholmblockid = stockholmblock.splitlines()[2].split()[-1]
            massive = True
        stockholmblockfh = StringIO(stockholmblock)
        outfasta = os.path.join(outdir, stockholmblockid + '.aln')
        records = SeqIO.parse(stockholmblockfh, "stockholm")
        print(records)
        if massive:
            fasta_fd = StringIO()
            SeqIO.write(records, fasta_fd, "fasta")
            fasta_fd.seek(0)
            origfasta_str = fasta_fd.read()
            parsedfasta_str = '>'.join(list(map(lambda x: '|'.join(x.split('|')[1:]), origfasta_str.split('>'))))
            with open(outfasta, 'w') as fastaout_fd:
                fastaout_fd.write(parsedfasta_str)
        else:
            SeqIO.write(records, outfasta, "fasta")
        seqtools.msa2fasta(outfasta, os.path.join(outdir, stockholmblockid + '.fa'))
        count += 1
    return count

def hmmpress_cleanup(globpath):
    globresults = glob.glob(globpath)
    if len(globresults):
        for old_hmmpress in glob.glob(globpath):
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


def rdrp_covfilter(domtblout, coverage):
    dfhmm = pd.read_table(domtblout, comment='#', header=None, delim_whitespace=True)
    dfhmm.columns = 'target_name target_accession tlen query_name query_accession qlen seq_evalue seq_score seq_bias domain_# domain_of domain_c-evalue domain_i-evalue domain_score domain_bias hmm_start hmm_end aln_start aln_end env_start env_end acc desc'.split()
    resultlist = []
    for (target, query), dftarget in dfhmm.groupby(['target_name', 'query_name']):
        covered = set()
        for domainid, record in dftarget.set_index('domain_#').iterrows():
            covered = covered | set(range(record['hmm_start'], record['hmm_end']+1))
            #print(len(covered))
        #print('result',len(covered))
        resultlist.append([target, query, record['seq_score'], len(covered), record['qlen'], len(covered)/record['qlen']])

    dfresultlist = pd.DataFrame(resultlist, columns=['target_name', 'query_name', 'seqscore','covered_length', 'hmm_length', 'coverage'])
    dfresultlist_rmdup = dfresultlist.sort_values('seqscore', ascending=False).drop_duplicates('target_name', keep='first')
    df_filtered = dfresultlist_rmdup[dfresultlist_rmdup.coverage >= coverage]
    return df_filtered


def seqkit_grep(seqid_list, fastain, fastaout):
    tmpid = str(int(time.time()*1e7))
    with open('.tmpseqid_%s.list' % tmpid, 'w') as fdout:
        fdout.writelines(list(map(lambda x: x+'\n',seqid_list)))
    (code, _, stderr) = runcommand("cat %s | sed 's/\/.*$//g' | seqkit grep -f .tmpseqid_%s.list > %s" % (fastain, tmpid, fastaout))
    if code:
        sys.stderr.write(curtime()+'[ERRO] seqkit exit with code %d\n' % code)
        sys.stderr.write(stderr.decode())
        sys.exit(1)
    else:
        os.remove('.tmpseqid_%s.list' % tmpid)


def preprocess(query, output, hmmdb):
    builded_hmm_path = os.path.join(output, 'builded_hmm')
    if not os.path.exists(builded_hmm_path):
        os.mkdir(builded_hmm_path)
    with open(hmmdb) as fdin:
        with open(os.path.join(builded_hmm_path, 'iter0_full.hmm'), 'w') as fdout:
            for line in fdin.readlines():
                if line.startswith('NAME') and not line.strip().endswith('_iter0'):
                    fdout.write(line.strip() + '_iter0\n')
                else:
                    fdout.write(line)
                
    clean_query = os.path.join(output, query.split('/')[-1])
    seqtools.clean_fasta(query, clean_query)
    sys.stderr.write(curtime()+'[INFO] sequence header cleanup done\n')
    return clean_query
    

def iterative_hmmsearch(query, wdir, ncpu, iteration_num, evalue, motifcov):
    sys.stderr.write('-----------------------Iterative_hmmsearch-------------------------\n')

    for i in range(1, iteration_num+1):
        sys.stderr.write('##Iteration %d\n' % i)
        iter_wdir = os.path.join(wdir, 'Iteration_%d' % i)
        if not os.path.exists(iter_wdir):
            os.mkdir(iter_wdir)
        else:
            sys.stderr.write(curtime()+'[WARN] Path Iteration_%d/ already exists! Files will be overwritten!\n' % i)

        #Step00. hmmpress
        _t0 = time.time()
        hmmpress_cleanup("{wdir}/builded_hmm/iter{last_iter}_full.hmm.*".format(wdir=wdir, last_iter=str(i-1)))
        code = os.system('hmmpress {wdir}/builded_hmm/iter{last_iter}_full.hmm > /dev/null 2>&1'.format(wdir=wdir, last_iter=str(i-1)))
        if code:
            sys.stderr.write(curtime()+'[ERRO] hmmpress exit with code %d\n'%code)
            sys.exit(1)
        _t1 = time.time()
        sys.stderr.write(curtime()+'[INFO] hmmpress done (Time Elapsed: %.3fs)\n' % (_t1-_t0))
        sys.stderr.write(curtime()+'[INFO] HMM search iteration %d running...\n' % i)

        #Step01. hmmsearch_run
        _t2 = time.time()
        code = os.system('hmmsearch --cpu {ncpu} -E {evalue} -o /dev/null -A {wdir}/Iteration_{iter}/hmmsearch.aln --tblout {wdir}/Iteration_{iter}/hmmsearch.tblout --domtblout {wdir}/Iteration_{iter}/hmmsearch.domtblout {wdir}/builded_hmm/iter{last_iter}_full.hmm {query}'.format(ncpu=ncpu, evalue=evalue, wdir=wdir, query=query, iter=i, last_iter=i-1))
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
        filtered_df = rdrp_covfilter(iter_wdir+'/hmmsearch.domtblout', motifcov)
        filtered_df.to_csv(os.path.join(iter_wdir,'profile_covstat.tsv'), sep='\t', index=None)
        
        #print(filtered_df)
        num_filtered = len(filtered_df['target_name'].drop_duplicates())
        _t5 = time.time()
        sys.stderr.write(curtime()+"[INFO] %d sequences covered >%d%% RdRp motif region (Time Elapsed: %.3fs)\n" % (num_filtered, int(100*motifcov), _t5-_t4))

        ##Step02.3 Fetch Confirmed Sequences
        num_profile, num_seq = 0, 0
        for profile, profile_df in filtered_df.groupby('query_name'):
            num_profile += 1
            seqlist = list(profile_df['target_name'].drop_duplicates())
            num_seq += len(seqlist)

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
        for file in bar(glob.glob(iter_wdir + '/sto2fas/*.filt.fa')):
            fileprefix = '.'.join(file.split('/')[-1].split('.')[:-2])
            code = os.system('cd-hit -d 0 -c 0.8 -aS 0.5 -i {query} -o {result_prefix} -T 4 > {result_prefix}.cdhit.log 2>&1'.format(query=file, result_prefix=os.path.join(cdhit_path, fileprefix)))
            if code:
                sys.stderr.write(curtime()+'[ERRO] cd-hit exit with code %s\n' % code)
                sys.exit(1)
        _t7 = time.time()
        sys.stderr.write(curtime()+'[INFO] cd-hit done (Time Elapsed: %.3fs)\n' % (_t7-_t6))

        ##Step03.2 Fetch Clustered Sequences
        hmmbuild_msa_path = os.path.join(iter_wdir, 'hmmbuild_MSA')
        try:
            os.mkdir(hmmbuild_msa_path)
        except FileExistsError:
            sys.stderr.write(curtime()+'[WARN] path %s exists, content will be overwritten!\n' % hmmbuild_msa_path)
        _num_seq = 0
        for _count, file in enumerate(glob.glob(cdhit_path + '/*.clstr')):
            dfclus = seqtools.read_cdhit(file, '0.8')
            profile = '.'.join(file.split('/')[-1].split('.')[:-1])
            keepseqid = list(dfclus.loc[dfclus['Cluster_Sign_0.8']=='*','Sequence_ID'])
            _num_seq += len(keepseqid)

            seqkit_grep(keepseqid, os.path.join(iter_wdir,'sto2fas', profile+'.aln'), os.path.join(hmmbuild_msa_path, '_'.join(profile.split('_')[:-1])+'_iter%d.aln'%i))

        _t8 = time.time()
        sys.stderr.write(curtime()+'[INFO] %d MSA seqs for %d profiles fetched (Time Elapsed: %.3fs)\n' % (_num_seq, _count+1, _t8-_t7))

        ##Step04. hmmbuild
        hmmbuild_path = os.path.join(iter_wdir, 'hmmbuild')
        try:
            os.mkdir(hmmbuild_path)
        except FileExistsError:
            sys.stderr.write(curtime()+'[WARN] path %s exists, content will be overwritten!\n' % hmmbuild_path)
        bar2 = progressbar.ProgressBar(widgets=[curtime(), '[INFO] hmmbuild running on profile ',progressbar.SimpleProgress(), ' (', progressbar.Percentage(), ')'])
        for file in bar2(glob.glob(hmmbuild_msa_path+'/*')):
            fileprefix = file.split('/')[-1]
            code = os.system('hmmbuild {hmm_prefix}.hmm {msa_prefix} > {hmm_prefix}.log 2>&1'.format(hmm_prefix=os.path.join(hmmbuild_path, fileprefix), msa_prefix=os.path.join(hmmbuild_msa_path, fileprefix)))
            if code:
                sys.stderr.write(curtime()+'[ERRO] hmmbuild exit with code %s\n' % code)
                sys.exit(1)

        code = os.system('cat {hmmbuild_path}/*.hmm > {iterative_hmm}'.format(hmmbuild_path=hmmbuild_path, iterative_hmm=os.path.join(wdir, 'builded_hmm', 'iter%d_full.hmm'%i)))
        sys.stderr.write(curtime()+'[INFO] Full iterative hmm file written to %s\n' % (os.path.join(wdir, 'builded_hmm', 'iter%d_full.hmm'%i)))
        _t9 = time.time()
        sys.stderr.write(curtime()+'[INFO] hmmbuild done (Time Elapsed: %.3fs)\n' % (_t9-_t8))
        sys.stderr.write(curtime()+'[INFO] Iteration %d finished (Time Elapsed: %.3fs)\n\n\n' % (i, _t9-_t1))
        

if __name__ == '__main__':
    main()

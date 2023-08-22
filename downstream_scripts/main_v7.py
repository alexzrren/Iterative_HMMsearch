#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 10:29:15 2023

@author: zirui
"""
from multiprocessing.pool import ThreadPool as Pool
import pandas as pd
from tqdm import tqdm
import numpy as np
import time
import re
import os


# %%

def process_vhdb():
    bx_vhdb_cols = 'qseqid qlen vhdb_sseqid vhdb_slen vhdb_qstart vhdb_qend vhdb_sstart vhdb_send vhdb_evalue vhdb_bitscore vhdb_length vhdb_pident vhdb_mismatch vhdb_gaps vhdb_stitle vhdb_qcovhsp vhdb_scovhsp'.split()
    blastx_data = pd.read_table('downstream_full/blast/vhdb/iterviralcontigs_vhdbcds.blastx', header=None, names=bx_vhdb_cols)

    blastx_data_rmdup = blastx_data.sort_values('vhdb_bitscore', ascending=False).drop_duplicates('qseqid', keep='first')
    blastx_data_rmdup['vhdb_svname'] = blastx_data_rmdup.vhdb_stitle.apply(lambda x: ' '.join(x.split('|')[0].split()[1:]))
    blastx_data_rmdup['vhdb_sprotein'] = blastx_data_rmdup.vhdb_stitle.apply(lambda x: x.split('|')[1])
    blastx_data_rmdup['vhdb_shost'] = blastx_data_rmdup.vhdb_stitle.apply(lambda x: x.split('|')[2])

    def find_viraltaxonomy(suffix, string):
        try:
            return re.search(suffix, string).group()
        except:
            return '-'

    def find_ifShiM(string):
        if re.search('ShiM', string):
            return True
        else:
            return False

    blastx_data_rmdup['vhdb_svorder'] = blastx_data_rmdup.vhdb_stitle.apply(lambda x: find_viraltaxonomy(r'\w+virales', x.split('|')[3]))
    blastx_data_rmdup['vhdb_svfamily'] = blastx_data_rmdup.vhdb_stitle.apply(lambda x: find_viraltaxonomy(r'\w+viridae', x.split('|')[3]))
    blastx_data_rmdup['vhdb_sifShiMang'] = blastx_data_rmdup.vhdb_stitle.apply(lambda x: find_ifShiM(x.split('|')[3]))
    blastx_data_rmdup['vhdb_vertebrate'] = blastx_data_rmdup.vhdb_stitle.apply(lambda x: 'Vertebrata' in x)

    vhdb_concise = blastx_data_rmdup[['qseqid','vhdb_svorder','vhdb_svfamily','vhdb_vertebrate']]
    return blastx_data_rmdup, vhdb_concise

# %% BLASTX NR

def process_nr():
    nr_vertebrate_vflist = '''Hepeviridae
Togaviridae
Picobirnaviridae
Arenaviridae
Hantaviridae
Nairoviridae
Peribunyaviridae
Phenuiviridae
Flaviviridae
Amnoonviridae
Orthomyxoviridae
Nodaviridae
Bornaviridae
Filoviridae
Paramyxoviridae
Pneumoviridae
Rhabdoviridae
Arteriviridae
Coronaviridae
Tobaniviridae
Caliciviridae
Picornaviridae
Sedoreoviridae
Spinareoviridae
Astroviridae
Adenoviridae
Asfarviridae
Anelloviridae
Circoviridae
Herpesviridae
Parvoviridae
Poxviridae
Papillomaviridae
Polyomaviridae'''.splitlines()

    blastnr_cols = ['qseqid', 'qlen'] + list(map(lambda x: 'nr_'+x, 'sseqid slen qstart qend sstart send evalue bitscore length pident mismatch gaps staxids sscinames sskingdoms skingdoms sphylums stitle qcovhsp scovhsp'.split()))
    blastnr_data = pd.read_table('downstream_full/blast/nr/iterviralcontigs_nr.blastx', header=None, names=blastnr_cols)

    blastnr_best75_datalist = []
    blastnr_best40_datalist = []
    blastnr_bestresult_datalist = []
    blastnr_fp_datalist = []
    for qseqid, qseqdf in tqdm(blastnr_data.groupby('qseqid')):
        #选择最高的比对分数50%/75%以内的结果
        qseqdf_best75 = qseqdf[qseqdf.nr_bitscore > qseqdf.nr_bitscore.max()*0.75]
        qseqdf_best40 = qseqdf[qseqdf.nr_bitscore > qseqdf.nr_bitscore.max()*0.40]
        blastnr_best75_datalist.append(qseqdf_best75)
        blastnr_best40_datalist.append(qseqdf_best40)
        #最优比对（最高bitscore）的一个或多个结果
        qseqdf_bestscore = qseqdf[qseqdf.nr_bitscore == qseqdf.nr_bitscore.max()]
        qseqdf_bestscore = qseqdf_bestscore.sort_values(['nr_pident', 'nr_scovhsp', 'nr_mismatch'], ascending=[False, False, True])
        qseqdf_bestscore_dropnil = qseqdf_bestscore[qseqdf_bestscore.nr_sskingdoms != '0']  #去除未知分类的结果
        if not len(qseqdf_bestscore_dropnil):
            blastnr_bestresult_datalist.append(qseqdf_bestscore.reset_index().iloc[0,:]) #如果没有非0分类的结果，就用最优比对结果
            continue
        else:
            blastnr_fp = qseqdf_bestscore_dropnil[(qseqdf_bestscore_dropnil.dropna().nr_sskingdoms.apply(lambda x: 'Viruses' not in x)) & (qseqdf_bestscore_dropnil.nr_length >= 200)]
            blastnr_bestresult_datalist.append(qseqdf_bestscore_dropnil.reset_index().iloc[0,:])
        blastnr_fp_datalist.append(blastnr_fp)

    blastnr_fp_full = pd.concat(blastnr_fp_datalist).sort_values('nr_bitscore', ascending=False).drop_duplicates('qseqid', keep='first') #假阳性的结果
    blastnr_data_rmdup = pd.concat(blastnr_bestresult_datalist, axis=1).T.drop('index', axis=1).reset_index(drop=True)
    blastnr_data_rmdup['pair'] = blastnr_data_rmdup.qseqid + '=' + blastnr_data_rmdup.nr_bitscore.astype(str)

    blastnr_best40data = pd.concat(blastnr_best40_datalist)
    blastnr_best75data = pd.concat(blastnr_best75_datalist)

    blastnr_best40data['nr_staxids_clean'] = blastnr_best40data.dropna()['nr_staxids'].apply(lambda x: str(list(filter(lambda y: y>0, list(map(int, x.split(';')))))[0]))
    blastnr_best75data['nr_staxids_clean'] = blastnr_best75data.dropna()['nr_staxids'].apply(lambda x: str(list(filter(lambda y: y>0, list(map(int, x.split(';')))))[0]))
    with open('.tmpfile.taxid', 'w') as fd:
        fd.writelines(list(blastnr_best40data['nr_staxids_clean'].dropna().drop_duplicates() +'\n'))

    os.system(r'cat .tmpfile.taxid | taxonkit lineage | taxonkit reformat -f "{o}\t{f}\t{g}" > .taxonkit.result.tsv')
    nr_vtaxinfo = pd.read_table('.taxonkit.result.tsv', header=None, names='nr_staxids_clean nr_slineage nr_svorder nr_svfamily nr_svgenus'.split(), sep='\t')
    os.remove('.tmpfile.taxid')
    os.remove('.taxonkit.result.tsv')

    nr_vtaxinfo['nr_staxids_clean'] = nr_vtaxinfo['nr_staxids_clean'].astype(str)
    blastnr_addtaxinfo40 = pd.merge(blastnr_best40data, nr_vtaxinfo, on='nr_staxids_clean', how='left')
    blastnr_addtaxinfo75 = pd.merge(blastnr_best75data, nr_vtaxinfo, on='nr_staxids_clean', how='left')

    def votevf(input_tuple):
        qseqid, qseqdf = input_tuple
        try:
            best_classified_bitscore = qseqdf.dropna(subset=['nr_svfamily']).reset_index().loc[0, 'nr_bitscore']  #获取最高比对分数的有分类结果的比对分数
        except:
            return [qseqid, '-', '-'] #如果在范围内都没有分类结果，则返回空并continue
        votedf = qseqdf[(qseqdf.nr_bitscore <= best_classified_bitscore) & (qseqdf.nr_bitscore > best_classified_bitscore*0.75)].dropna(subset=['nr_svfamily']) #选择最高比对分数的75%以内的有分类结果
        bitscore_voting = votedf.groupby('nr_svfamily')['nr_bitscore'].sum().sort_values(ascending=False) #对分类结果进行投票
        bitscore_sum = votedf['nr_bitscore'].sum() #计算总比对分数
        return [qseqid, bitscore_voting.index[0], bitscore_voting[bitscore_voting.index[0]]/bitscore_sum] #返回投票结果(即投票获得的bitscore最高的)

    time1 = time.time()
    input_tuplelist = blastnr_addtaxinfo40.groupby('qseqid')
    pool = Pool(processes=16)
    result = pool.map(votevf, input_tuplelist)
    pool.close()
    pool.join()
    time2 = time.time()

    votingdf = pd.DataFrame(result, columns=['qseqid', 'nr_vote_vfamily', 'approval_rate'])

    blastnr_filtedtaxinfo = blastnr_addtaxinfo75.dropna(axis=0, how='any', subset=['nr_slineage']).dropna(axis=0, how='all', subset=['nr_svorder', 'nr_svfamily'])
    nr_75svorderlist = blastnr_filtedtaxinfo.groupby('qseqid').apply(lambda x: ';'.join(set(x.nr_svorder.dropna())))
    nr_75svordercount = blastnr_filtedtaxinfo.groupby('qseqid').apply(lambda x: len(set(x.nr_svorder.dropna())))
    nr_75svfamilylist = blastnr_filtedtaxinfo.groupby('qseqid').apply(lambda x: ';'.join(set(x.nr_svfamily.dropna())))
    nr_75svfamilycount = blastnr_filtedtaxinfo.groupby('qseqid').apply(lambda x: len(set(x.nr_svfamily.dropna())))
    nr_taxstat = pd.concat([nr_75svordercount, nr_75svorderlist, nr_75svfamilycount, nr_75svfamilylist], axis=1)
    nr_taxstat.columns = ['nr_75svordercount', 'nr_75svorderlist', 'nr_75svfamilycount', 'nr_75svfamilylist']

    blastnr_filtedtaxinfo_rmdup = blastnr_filtedtaxinfo.sort_values('nr_bitscore', ascending=False).drop_duplicates('qseqid', keep='first')

    blastnr_filtedtaxinfo_rmdup['pair'] = blastnr_filtedtaxinfo_rmdup.qseqid + '=' + blastnr_filtedtaxinfo_rmdup.nr_bitscore.astype(str)
    blastnr_filtedtaxinfo_rmdup['nr_validtaxinfo'] = True

    lack = set(blastnr_data_rmdup.qseqid) - set(blastnr_filtedtaxinfo_rmdup.qseqid)
    blastnr_full = pd.concat([blastnr_filtedtaxinfo_rmdup, blastnr_data_rmdup[blastnr_data_rmdup.qseqid.isin(lack)]])
    blastnr_full['nr_validtaxinfo'].fillna(False, inplace=True)
    blastnr_full['nr_besthit'] = False
    blastnr_full.loc[blastnr_full.pair.isin(blastnr_data_rmdup.pair), 'nr_besthit'] = True
    blastnr_full = blastnr_full.drop(['pair'], axis=1)
    blastnr_full['nr_FP'] = False
    blastnr_full.loc[blastnr_full.qseqid.isin(blastnr_fp_full.qseqid), 'nr_FP'] = True

    blastnr_full = pd.merge(blastnr_full, votingdf, on='qseqid', how='left')
    blastnr_full = pd.merge(blastnr_full, nr_taxstat.reset_index(), on='qseqid', how='left')

    blastnr_full['nr_vfamily_add75'] = blastnr_full['nr_svfamily']
    blastnr_full.loc[(blastnr_full['nr_75svfamilycount']>0) & (blastnr_full['nr_vfamily_add75'].isna()) & (blastnr_full['nr_75svfamilycount']<2), 'nr_vfamily_add75'] = blastnr_full['nr_75svfamilylist']

    blastnr_full['nr_vertebrate'] = blastnr_full['nr_vote_vfamily'].isin(nr_vertebrate_vflist)

    nr_concise = blastnr_full[['qseqid', 'nr_FP', 'nr_sseqid' ,'nr_svorder', 'nr_svfamily', 'nr_svgenus', 'nr_sscinames', 'nr_besthit', 'nr_pident', 'nr_vfamily_add75', 'nr_vote_vfamily' , 'approval_rate', 'nr_vertebrate']].fillna('-')
    return blastnr_full, nr_concise



# %% BLASTN NT

def process_nt():
    blastnt_cols = ['qseqid', 'qlen'] + list(map(lambda x: 'nt_'+x, 'sseqid slen qstart qend sstart send evalue bitscore length pident mismatch gaps staxid ssciname scomname sblastname sskingdom staxids sscinames scomnames sblastnames sskingdoms stitle qcovs qcovhsp qcovus'.split()))
    blastnt_data = pd.read_table('downstream_full/blast/nt/iterviralcontigs_nt.blastn', header=None, names=blastnt_cols)
    blastnt_data_rmdup = blastnt_data.sort_values('nt_bitscore', ascending=False).drop_duplicates('qseqid', keep='first')
    blastnt_fp = blastnt_data_rmdup[(blastnt_data_rmdup.nt_sskingdoms!='Viruses') & (blastnt_data_rmdup.nt_length >= 200)]
    blastnt_data_rmdup['nt_FP'] = False
    blastnt_data_rmdup.loc[blastnt_fp.index, 'nt_FP'] = True

    nt_concise = blastnt_data_rmdup[['qseqid', 'nt_FP']]
    return blastnt_data_rmdup, nt_concise

# %% motif coverage filtering

def filt_motifcov(itersearch_path):
    hmmsearch_cols_clean = 'target_name target_accession tlen query_name query_accession qlen seq_evalue seq_score seq_bias domain_# domain_of domain_c-evalue domain_i-evalue domain_score domain_bias hmm_start hmm_end aln_start aln_end env_start env_end acc desc'.split()

    iterhmmsearchlist = []
    for i in range(1,11):
        hmmsearch_result = pd.read_table(os.path.join(itersearch_path, 'Iteration_%d/hmmsearch.domtblout'%i), comment='#', delim_whitespace=True, names=hmmsearch_cols_clean)
        iterhmmsearchlist.append(hmmsearch_result)

    fulliter = pd.concat(iterhmmsearchlist)
    fulliter['iteration'] = fulliter['query_name'].apply(lambda x: x[-1]).astype(int) + 1
    fulliter['target_name_withquery'] = fulliter['target_name'] + '_' + fulliter['query_name']

    covresult = []
    keep_bestlist = fulliter.sort_values('seq_score', ascending=False).drop_duplicates('target_name')['target_name_withquery']
    fulliter_bestresult = fulliter[fulliter['target_name_withquery'].isin(keep_bestlist)]
    for seq, seqdf in fulliter_bestresult.groupby('target_name'):
        profilemin = 99999999999
        profilemax = 0
        hmmlen = seqdf.loc[seqdf.index[0],'qlen']
        seqvec = np.zeros(hmmlen)
        for index, record in seqdf.iterrows():
            seqvec[record['hmm_start']:record['hmm_end']] = 1
            if record['hmm_start'] < profilemin:
                profilemin = record['hmm_start']
            if record['hmm_end'] > profilemax:
                profilemax = record['hmm_end']
        cov = np.count_nonzero(seqvec) / seqvec.size
        covresult.append([seq, record['query_name'], profilemin, profilemax, cov])


    fulliter_domaincov_sum = pd.DataFrame(covresult, columns=['qorfid', 'besthit_hmm', 'Profile_start', 'Profile_end', 'Profile_cov'])
    fulliter_domaincov_sum['qseqid'] = fulliter_domaincov_sum['qorfid'].apply(lambda x: '_'.join(x.split('_')[:-1]))
    fulliter_domaincov_sum = fulliter_domaincov_sum.sort_values('Profile_cov', ascending=False).drop_duplicates('qseqid')
    fulliter_domaincov_sum['valid_Profile_cov'] = False
    fulliter_domaincov_sum.loc[fulliter_domaincov_sum.Profile_cov >= 0.50, 'valid_Profile_cov'] = True
    #fulliter_domaincov_sum = fulliter_domaincov_sum[['qseqid', 'target_name', 'Profile_cov', 'valid_Profile_cov']]
    bestorflist = list(fulliter_domaincov_sum['qorfid'])
    return fulliter_domaincov_sum, bestorflist

# %% hmmsearch RdRp-scan

def hmm_classify(bestorflist):
    hmmsearch_cols = list(map(lambda x: 'rdrpscan_'+x, 'target_name target_accession tlen query_name query_accession qlen seq_evalue seq_score seq_bias domain_# domain_of domain_c-evalue domain_i-evalue domain_score domain_bias hmm_start hmm_end aln_start aln_end env_start env_end acc desc'.split()))
    rdrpscan = pd.read_table('downstream_full/hmmsearch_RdRpscan/hmmsearch_RdRpscan.domtblout', comment='#', delim_whitespace=True, names=hmmsearch_cols)
    rdrpscan = rdrpscan[rdrpscan['rdrpscan_target_name'].isin(bestorflist)]
    rdrpscan['qseqid'] = rdrpscan['rdrpscan_target_name'].apply(lambda x: '_'.join(x.split('_')[:-1]))
    rdrpscan['rdrpscan_svorder'] = rdrpscan['rdrpscan_query_name'].apply(lambda x: x.split('.')[1].split('_')[0])
    rdrpscan_besthit = rdrpscan.sort_values('rdrpscan_seq_score', ascending=False).drop_duplicates('rdrpscan_target_name')
    #stat = rdrpscan_besthit.value_counts('qseqid')

    rdrpscan_concise = rdrpscan_besthit[['qseqid', 'rdrpscan_svorder', 'rdrpscan_seq_evalue']]
    rdrpscan_full = rdrpscan_besthit[['qseqid', 'rdrpscan_svorder', 'rdrpscan_seq_evalue', 'rdrpscan_seq_score']]
    return rdrpscan_full, rdrpscan_concise


# %% RdRp AAI90 cluster

def read_cdhit(file, identity, length=False):
    clusinfo = []
    with open(file) as clus:
        for line in clus.readlines():
            if line.startswith('>'):
                clusterid = line.strip().split(' ')[-1]
            else:
                if line != '':
                    seqid = line.split(' ')[1].lstrip('>').rstrip('...')
                    repsign = line.split(' ')[-1].strip()
                    seqlen = line.split('\t')[1].split(',')[0].rstrip('nt')
                    clusinfo.append([clusterid, seqid, seqlen, repsign])
    identity = str(identity)
    dfclus = pd.DataFrame(clusinfo, columns=['Cluster_AAI'+identity, 'qorfid', 'Sequence_Length','ClusterSign_AAI'+identity])
    if length:
        return dfclus
    else:
        return dfclus.drop('Sequence_Length', axis=1)


def process_clusinfo(bestorflist):
    rcr90 = read_cdhit('downstream_full/sequences/putative_markerproteins/AAI90.clstr', 90, length=False)
    rcr90['qseqid'] = rcr90['qorfid'].apply(lambda x: '_'.join(x.split('_')[:-1]))
    rcr90 = rcr90[rcr90['qorfid'].isin(bestorflist)]

    rcr90_stat = rcr90.value_counts('Cluster_AAI90')
    rcr90_stat = pd.DataFrame(rcr90_stat, columns = ['AAI90_Seqcount']).reset_index()
    rcr90_stat.columns = ['Cluster_AAI90', 'AAI90_Seqcount']
    rcr90 = pd.merge(rcr90, rcr90_stat, on='Cluster_AAI90', how='left')

    clus30 = pd.read_table('downstream_full/sequences/putative_markerproteins/ident30cov60_clus_cluster.tsv', header=None, names=['30aai_rep', 'qorfid'])
    clus30 = clus30[clus30['qorfid'].isin(bestorflist)]
    clus30_id2rep = pd.DataFrame([[str(i),rep] for i,rep in enumerate(set(clus30['30aai_rep']))], columns=['Cluster_AAI30', '30aai_rep'])
    clus30 = pd.merge(clus30_id2rep, clus30, on='30aai_rep', how='left')
    clus30['ClusterSign_AAI30'] = '-'
    clus30.loc[clus30['30aai_rep']==clus30['qorfid'], 'ClusterSign_AAI30'] = '*'
    clus30['qseqid'] = clus30['qorfid'].apply(lambda x: '_'.join(x.split('_')[:-1]))

    clusresult = pd.merge(rcr90[['qseqid', 'Cluster_AAI90','ClusterSign_AAI90', 'AAI90_Seqcount']], clus30[['qseqid', 'Cluster_AAI30', 'ClusterSign_AAI30']], on='qseqid', how='left')
    return clusresult


# %% summaryout

if __name__=='__main__':
    seqlen = pd.read_table('downstream_full/sequences/putative_viralcontigs/putative_viralcontigs.lenstat.tsv', header=None, names=['qseqid', 'qlen'])
    seqlen['viraltype'] = 'RNA'
    with open('downstream_full/sequences/putative_viralcontigs/dna.list') as fd:
        dna_contigs = fd.read().splitlines()
    seqlen.loc[seqlen['qseqid'].isin(dna_contigs), 'viraltype'] = 'DNA'

    nr_full, nr_concise = process_nr()
    nt_full, nt_concise = process_nt()
    vhdb_full, vhdb_concise = process_vhdb()
    rna_domcov, rna_bestorflist = pd.DataFrame(), []
    dna_domcov, dna_bestorflist = filt_motifcov('output_dna_v1')
    bestorflist = rna_bestorflist + dna_bestorflist
    fulliter_domaincov_sum = pd.concat([rna_domcov, dna_domcov], axis=0)
#    rdrpscan_full, rdrpscan_concise = hmm_classify(bestorflist)
#    clusresult = process_clusinfo(bestorflist)
    summary_concise0 = pd.merge(seqlen, fulliter_domaincov_sum, on='qseqid', how='left')
    summary_concise2 = summary_concise0
#    summary_concise1 = pd.merge(summary_concise0, clusresult, on='qseqid', how='left')
#    summary_concise2 = pd.merge(summary_concise1, rdrpscan_concise, on='qseqid', how='left')
    summary_concise3 = pd.merge(summary_concise2, nt_concise, on='qseqid', how='left')
    summary_concise4 = pd.merge(summary_concise3, nr_concise, on='qseqid', how='left')
    summary_concise5 = pd.merge(summary_concise4, vhdb_concise, on='qseqid', how='left')
    summary_concise5.to_csv('Marine_DNAvirus.csv', index=False)
    summary_clean = summary_concise5[(summary_concise5.valid_Profile_cov) & (~summary_concise5.nt_FP.fillna(False)) & (~summary_concise5.nr_FP.fillna(False)) ]
    summary_clean_vert = summary_clean[(summary_clean.nr_vertebrate) & (summary_clean.vhdb_vertebrate.fillna(True))]
    summary_clean_vert.to_csv('Marine_TPDNA_Vertvirus', index=False)
    
    summary_full2 = summary_concise0
    summary_full3 = pd.merge(summary_full2, nt_full.drop(['qlen','nt_qcovs', 'nt_qcovhsp', 'nt_qcovus'], axis=1), on='qseqid', how='left')
    summary_full4 = pd.merge(summary_full3, nr_full.drop('qlen', axis=1), on='qseqid', how='left')
    summary_full5 = pd.merge(summary_full4, vhdb_full.drop('qlen', axis=1), on='qseqid', how='left')
    summary_full5.to_csv('summary_full.csv', index=False)
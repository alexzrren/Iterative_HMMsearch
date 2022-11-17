def dict2fasta(fastadict, file):
    fasta = ''
    for seqid, seq in fastadict.items():
        fasta += '>' + seqid +'\n' 
        while len(seq) >= 60:
            fasta += seq[0:60] + '\n'
            seq = seq[60:]
        if seq == '':
            continue
        else:
            fasta += seq +'\n'
    with open(file, 'w') as fastaout:
        fastaout.write(fasta)
            

def msa2fasta(filein, fileout):
    msadict = fasta2dict(filein)
    fastadict = dict( (k,v.replace('-','')) for k,v in msadict.items() )
    dict2fasta(fastadict, fileout)  


def fasta2dict(file):
    with open(file) as fasta:
        fastadict = dict()
        seqid = ''
        seq = ''
        for line in fasta.read().splitlines():
            if line.startswith('>'):
                fastadict[seqid] = seq
                seq = ''
                seqid = line.lstrip('>').split()[0].split('/')[0]
            else:
                seq += line
    fastadict[seqid] = seq
    del fastadict['']
    return fastadict


def clean_fasta(file, cleanfile):
    import os
    os.system('awk \''+'{'+'print $1'+'}'+'\' '+'{file} > {cleanfile}'.format(file=file, cleanfile=cleanfile))
    
    
def read_cdhit(file, identity, length=False):
    """
    This function is to read the result file (cdhit.clstr) of CD-HIT, and convert it in to pandas DataFrame
    :param file: CD-HIT聚类结果文件(.clstr)
    :param identity: 用于CD-HIT聚类的一致性
    :param length: 是否返回序列长度列
    :return: 返回值是存储CD-HIT聚类信息的pandas.DataFrame
    """
    import pandas as pd
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
    dfclus = pd.DataFrame(clusinfo, columns=['Cluster_'+identity, 'Sequence_ID', 'Sequence_Length','Cluster_Sign_'+identity])
    if length:
        return dfclus
    else:
        return dfclus.drop('Sequence_Length', axis=1)
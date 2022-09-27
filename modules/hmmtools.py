#!/home/renzirui/miniconda3/bin/python

import sys

def read_hmm(input = 'stdin'):
    if input == 'stdin':
        sys.stderr.write('[INFO] stdin mode on')
        content = sys.stdin
    else:
        with open(input, 'r') as fd:
            content = fd.readlines()
    buffer = ''
    contentblocks = []
    for line in content:
        if line.startswith('//'):
            buffer += line
            contentblocks.append(buffer)
            buffer = ''
        else:
            buffer += line
    sys.stderr.write('[INFO] Successfully read %s hmm profiles\n' % (len(contentblocks)))
    return contentblocks


def block2dict(contentblocks):
    contentdict = dict((i.splitlines()[1].split()[-1], i) for i in contentblocks)
    return contentdict


def fetch_profile(contentblocks, name):
    try:
        sys.stdout.write(block2dict(contentblocks)[name])
    except:
        sys.stderr.write('[ERRO] Profile name %s not found\n' % (name))


def hmm2seq_core(contentblock):
    matchsign = 1
    msa_info = []
    for line in contentblock.splitlines()[18:]:
        if line.split()[0] != str(matchsign):
            continue
        else:
            matchsign += 1
            msa_info.append(line.split()[-5:-1])

    return ''.join([ i[1] for i in msa_info ])
    

def cut(obj, sec):
    return [obj[i:i+sec] for i in range(0,len(obj),sec)]


def hmm2seq(contentblocks, output='return'):
    if output == 'stdout':
        stdout = True
    elif output == 'return':
        stdout = False
    else:
        print('[ERRO] illegal output format')

    contentdict = block2dict(contentblocks)
    if stdout:
        for name, contentblock in contentdict.items():
            seq = hmm2seq_core(contentblock)
            sys.stderr.write('[INFO] Profile %s converted\n' % (name))
            sys.stdout.write('>' + name + '\n')
            for seqline in cut(seq, 60):
                sys.stdout.write(seqline.upper() + '\n')
        return None
    else:
        return_seq = ''
        for name, contentblock in contentdict.items():
            seq = hmm2seq_core(contentblock)
            sys.stderr.write('[INFO] Profile %s converted\n' % (name))
            return_seq += ('>' + name + '\n')
            for seqline in cut(seq, 60):
                return_seq += (seqline.upper() + '\n')


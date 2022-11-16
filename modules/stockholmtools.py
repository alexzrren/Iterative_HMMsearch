from Bio import SeqIO
import sys

def read_stockholm():
    buffer = ''
    contentblocks = []
    for line in sys.stdin:
        if line.startswith('//'):
            buffer += line
            contentblocks.append(buffer)
            buffer = ''
        else:
            buffer += line
    sys.stderr.write('[INFO] Successfully read MSA of %d profiles\n' % (len(contentblocks)))
    return contentblocks


def block2dict(contentblocks):
    contentdict = dict((i.splitlines()[1].split()[-1],i) for i in contentblocks)
    return contentdict


def fetch_profile(contentblocks, name):
    try:
        sys.stdout.write(block2dict(contentblocks)[name])
    except:
        sys.stderr.write('[ERRO] Profile name %s not found\n' % (name))


def parse_profile(contentblock):
    profiledict = dict()
    for line in contentblock.splitlines():
        pass


def stockholm2fasta(stockholm, fasta):
    records = SeqIO.parse(stockholm, "stockholm")
    count = SeqIO.write(records, fasta, "fasta")
    print("Converted %i records" % count)

if __name__ == '__main__':
    contentblocks = read_stockholm()
    fetchname = sys.argv[1]
    #fetch_profile(contentblocks, fetchname)
    parse_profile(block2dict(contentblocks)[fetchname])
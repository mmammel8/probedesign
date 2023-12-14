import os
import sys
import re
import gzip
from subprocess import Popen, PIPE
from string import digits

outfile = ""
fnameIn = "cx.fasta.gz"
#list of accession numbers to include
flist = 'salm_list.txt'
serovar = "Salmonella_enterica_Newport_cluster3"
strain_dict = dict()
nvar = 0

#BLAST databases to use
db = 'SALM EXCL2 PLSM'

#number of genomes with <>1 matches allowed
MAXPOINTS = 5
PROBELEN = 120
MAXPROBE = 80000
MAXOUTS = 40

#minimum match percent
MINPERC = 90.0

#minimum match length
MINLEN = PROBELEN * 9 // 10

#number of processors for BLAST
NUM_THREADS = 8

#must include fasta file of all accessions used in fgenes file
#saved in same directory, and named accession.fna

#globals
blin = "seq1.txt"
blout = "result1.txt"

def capture(var):
    """
    Read BLAST result file
    """
    tallyin = 0
    tallyout = 0
    data_file = open(blout, 'r')
    data = ''.join(data_file.readlines())
    data_file.close()
    lines = data.split('\n')
    for line in lines:
        row = line.rstrip().split('\t')
        if len(row) > 11:
            query = row[0]
            subject = row[1]
            percent = float(row[2])
            alen = int(row[3])
            bscr = float(row[11])
            pos = subject.find("_")
            if pos > -1:
                 subject = subject[:pos]
            if subject in strain_dict:
                var2 = strain_dict[subject]
            else:
                var2 = "exclude"
            if percent >= MINPERC and alen >= MINLEN:
                if var2 == var:
                    tallyin += 1
                else:
                    tallyout += 1
    return (tallyin, tallyout)

def blast_it():
    """
    perform BLAST command
    """
    blast_args = ['blastn', '-task', 'megablast', '-db', db, '-query', blin, '-out', blout, '-max_target_seqs', '20000', '-outfmt' , '6', '-num_threads', str(NUM_THREADS)]
    process = Popen(blast_args, stdout=PIPE)
    (stdout, stderr) = process.communicate()


def rev_comp(seq):
    """
    return reverse complement
    """
    complement = {"A":"T", "C":"G", "G":"C", "T":"A"}
    return "".join([complement[base] for base in reversed(seq)])

def save_q(seq, locus):
    """
    save seq of locus to file for BLAST
    """
    data_file = open(blin, 'w')
    out = ">" + locus + "\n"
    data_file.write(out)
    out = seq + "\n"
    data_file.write(out)
    data_file.close()

def process(refseq, contig, var):
    """
    blast probe and record results
    """
    llen = len(refseq)
    res = []
    for start in range(0, llen - PROBELEN, PROBELEN):
        seq = refseq[start : start+PROBELEN]
        label = contig + "_" + str(start + 1)
        save_q(seq, label)
        blast_it()
        tallyin, tallyout = capture(var)
        score = 100000 - abs(tallyin - nvar) * 10 - tallyout*100
        rec = (score, label, var, nvar, tallyin, tallyout)
        res.append(rec)
    return res

if __name__ == '__main__':
    #get list of strains to include
    ref_acc = None
    data_file = open(flist, 'r')
    data = ''.join(data_file.readlines())
    data_file.close()
    lines = data.split('\n')
    for line in lines:
        row = line.split('\t')
        if len(row) > 1:
            strain_dict[row[0]] = row[1]
            if row[1] == serovar:
                #acc_set.add(row[0])
                nvar += 1
                if ref_acc == None:
                    ref_acc = row[0]
                    print("ref ", row[0])

    print(nvar, " ingroups")
    outfile = serovar + ".txt"
    seq = ""
    contig = ""
    acc = None
    result = []
    with gzip.open(fnameIn,'rt') as f:
    #with open(fnameIn,'r') as f:
        for line in f:
            if len(line) > 1:
                line = line.strip()
                if line.startswith('>'):
                    #header
                    if acc == ref_acc and seq != "":
                        res = process(seq, contig, serovar)
                        print(serovar, contig, len(res))
                        result.extend(res)
                    contig = line[1:]
                    acc = contig
                    pos = contig.find("_")
                    if pos > -1:
                        acc = contig[:pos]
                    seq = ""
                else:
                    seq += line.upper()
    f.close()
                           

t_file = open(outfile, 'w')
result.sort(reverse = True)
n = min(len(result), MAXPROBE)
for i in range(n):
    out = "\t".join(map(str,result[i])) + "\n"
    t_file.write(out)
t_file.close()





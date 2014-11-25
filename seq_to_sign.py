#! /sw/bin/python
import sys
import getopt
from pyfasta import Fasta
import string
complementTable = string.maketrans('ACTGRYSWKMBDHVNactgryswkmbdhvn','TGACYRSWMKVHDBNtgacyrswmkvhdbn')

def makeRC(origSeq):
    seqList = list(origSeq)
    seqList.reverse()
    seqRC = string.join(seqList,'')
    seqRC = seqRC.translate(complementTable)
    return seqRC

def make_seq_files(fasta_file, genomic_regions, window, name, motif):
    infile = open(genomic_regions, 'r')
    outfilename = name + '.fasta'
    print 'making' 
    print outfilename
    outfile = open(outfilename, 'w')
    dic_regions ={}
    chr = 0
    begin = 2
    end = 3
    strand = 1
    while 1:
        line = infile.readline()
        if not line: break
        if line.startswith('#'): continue
        linesplit = line.split()
        chrN = linesplit[chr]
        start = (int(linesplit[begin])) - int(window)
        stop = (int(linesplit[end])) + int(window)
        newfeature = (chrN, start, stop, string.split(linesplit[strand], motif))
        if chrN not in dic_regions:
            dic_regions[chrN]=[]
        dic_regions[chrN].append(newfeature)
    genome = Fasta(fasta_file)
    for key in dic_regions:
        for i in dic_regions[key]:
            chr = str(i[0])
            region_1 = int(i[1])
            region_2 = int(i[2])
            if i[3][0] == '+':
                seq = genome[chr][region_1-1:region_2+1]
                outfile.write('%s\n%s\n'%('>'+chr+':'+str(region_1+1)+'-'+str(region_2), seq))
            else:
                seq = genome[chr][region_1-2:region_2]
                outfile.write('%s\n%s\n'%('>'+chr+':'+str(region_1)+'-'+str(region_2-1), makeRC(seq)))
    infile.close()
    outfile.close()
    return

def makedf(fastaout):
    outfile=open(fastaout.split('.fasta')[0] + '.txt', 'w')
    f = open(fastaout)
    lines = f.readlines()
    f.close()
    grid = [[0]*(len(lines[1])-4) for i in range(len(lines)/2)]
    count = -1
    for i in lines:
        if i.startswith('>'): continue
        else:
            count = count + 1
            for j in range(len(i)-4):
                grid[count][j] =  i[j] + i[j+1] + i[j+2] + i[j+3]
    for k in range(0,len(grid)):
        h = string.join(grid[k],'\t')
        outfile.write('%s\n'%(h))
        
    outfile.close()
    return

def main(argv):
    try:
        opts, args = getopt.getopt(argv, "hi:f:w:n:m:", ["help", "input=", "fasta=", "window=", "name=", "motif="])
    except getopt.GetoptError, err:
        print str(err)
        sys.exit(2)
    namefile = False
    fast = False
    win = False
    name = False
    for opt, arg in opts:
        if opt in ('-i', '--input'):
            namefile = arg
        elif opt in ('-f', '--fasta'): 
            fast = arg
        elif opt in ('-w', '--window'): 
            win = arg
        elif opt in ('-n', '--name'): 
            name = arg
        elif opt in ('-m', '--motif'): 
            motif = arg
        elif opt in ('-h', '--help'):
            print 'python ~/seq_to_sign/seq_to_sign.py -i ~/seq_to_sign/test_files/mast.srf.hg19.chr1.txt -f ~/seq_to_sign/test_files/hg19.chr1.fa -w 20 -n SRF -m 1\n-i processed mast file in bed format\n-f genomic fasta file\n-w distance to look from motif center\n-n prefix\n-m motif number from mast file'
            sys.exit()
    if namefile and fast and win and name and motif:
        make_seq_files(fast, namefile, win, name, motif)
        newfast = name + '.fasta'
        makedf(newfast)
if __name__ == "__main__":
    main(sys.argv[1:])




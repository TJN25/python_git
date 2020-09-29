#!/usr/bin/python
import sys
from Bio import SeqIO

help = '''
    sequence_cleaner.py v 0.1 (September 2020) is a script for {}.
    Usage:
         sequence_cleaner.py [options] <fastaFile> <outFile>
    
    Options:
        -h	Display this help
        -q  Supress messages
    Input
        -i	<fastaFile> the input
        -o	<outFile> the output
        -l  <min_length> the minimum sequence length
        -u  <max_length> the maximum sequence length

    
'''

def usage():
    print help

def rungetopts():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "i:o:qh", ["fastaFile", "outFile", "quiet", "help"])
    except getopt.GetoptError as err:
        print(err)
        usage()
        sys.exit(2)
    fastaFile = ""
    outFile = ""
    for o, a in opts:
            if o in ("-h", "--help"):
                usage()
                sys.exit()
            elif o in ("-i", "--fastaFile"):
                fastaFile = a
            elif o in ("-o", "--outFile"):
                outFile = a
            elif o in ("-l", "--minlength"):
                min_length = a
            elif o in ("-u", "--maxlength"):
                max_length = a
            else:
                assert False, "unhandled option"
    if outFile == "":
        print "-o <outFile> missing. For more help use -h"
        sys.exit(2)
    if fastaFile == "":
        print "-i <fastaFile> missing. For more help use -h"
        sys.exit(2)
    return(fastaFile, outFile, min_length, max_length)


def sequence_cleaner(fastaFile, outFile, min_length=50, max_length=500, por_n=100):
    for seq_record in SeqIO.parse(fastaFile, "fasta"):
        iter += 1
        sequence = str(seq_record.seq).upper()
        if (
                len(sequence) >= min_length
                and (float(sequence.count("N")) / float(len(sequence))) * 100 <= por_n
                and len(sequence) <= max_length
        ):
            seq_name = seq_record.id
            sequence = seq_record.seq
            outFile.write(f">{seq_name}\n{sequence}\n")


def main():

    inFile, outFile, min_length, max_length = rungetopts()
    sequence_cleaner(fastaFile=fastaFile, outFile=outFile, min_length=min_length, max_length=max_length)



if __name__ == "__main__":
    main()

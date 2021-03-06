#!/usr/bin/python
import sys
from Bio import SeqIO
import getopt

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
        -l  <min_length> the minimum sequence length (default = 50)
        -u  <max_length> the maximum sequence length (default = 500)

    
'''

def usage():
    print help

def rungetopts():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "i:o:l:u:qh", ["fastaFile", "outFile", "minlength", "maxlength", "quiet", "help"])
    except getopt.GetoptError as err:
        print(err)
        usage()
        sys.exit(2)
    fastaFile = ""
    outFile = ""
    min_length = 50
    max_length = 500
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


def sequence_cleaner(fastaFile, outFile, min_length=50, max_length=500):
    # print("running sequence cleaner using: min_length = %s, max_length = %s" % (min_length, max_length))
    outFile = open(outFile, 'a')
    for seq_record in SeqIO.parse(fastaFile, "fasta"):
      #  print(seq_record.id)
        sequence = str(seq_record.seq).upper()
        # print("%s %s" % (len(sequence), min_length))
        if len(sequence) > int(min_length) and len(sequence) <= int(max_length):
            seq_name = seq_record.id
            nucleotides = seq_record.seq
            print("writing %s" % seq_name)
            outFile.write(">%s\n%s\n" %(seq_name, nucleotides))
        else:
            seq_name = seq_record.id
            nucleotides = seq_record.seq
            print("%s failed filter.\n %s\n" %(seq_name, nucleotides))
def main():
    # print("running getopts")
    fastaFile, outFile, min_length, max_length = rungetopts()
    sequence_cleaner(fastaFile=fastaFile, outFile=outFile, min_length=min_length, max_length=max_length)



if __name__ == "__main__":
    main()

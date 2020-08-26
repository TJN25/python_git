#!/usr/bin/python
import sys
from Bio import SeqIO
from Bio.Alphabet import generic_dna
import getopt
import os



help = '''
    get_sRNA_sequences.py v 0.1 (August 2020) is a script for taking predicted sRNA regions and extracting the DNA sequence from 
    fasta files.
    Usage:
     get_sRNA_sequences.py <accession>
    
    Options:
        -h	Display this help
        -q  Supress messages
    Input
        -a	<accession> the genome accession to use for the input and output files

    
'''

def usage():
    print help

def rungetopts():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "a:qh", ["accession", "quiet", "help"])
    except getopt.GetoptError as err:
        print(err)
        usage()
        sys.exit(2)
    accession = ""
    for o, a in opts:
            if o in ("-h", "--help"):
                usage()
                sys.exit()
            elif o in ("-a", "--accession"):
                accession = a
            else:
                assert False, "unhandled option"
    if accession == "":
        print "-a <accession> missing. For more help use -h"
        sys.exit(2)
    return(accession)




#
# def makeoutputdirectory(write_path):
#     if os.path.isdir(write_path) == False:
#         try:
#             os.mkdir(write_path)
#         except OSError:
#             print ("Creation of the directory %s failed" % accession)
#             sys.exit(2)
#     directory = os.listdir(write_path)
#     if len(directory) != 0:
#         print "Examples of files in %s" % write_path
#         print directory[0:4]
#         query_user = raw_input("%s is not an empty directory. Continue anyway y/n (this may write over existing files): " % write_path)
#         if query_user == "y":
#             print "Using %s as directory" % write_path
#         else:
#             print "Exiting script"
#             sys.exit(2)
#
# def concatenateSequence(fastaFile):
#     my_seq = fastaFile[0].seq
#     i = 0
#     for seq in fastaFile:
#         if i == 0:
#             i += 1
#             continue
#         i += 1
#         my_seq = my_seq + seq.seq
#     return my_seq
#
# def writeSequences(inFile,my_seq,accession,write_path):
#     i = 0
#     for line in inFile:
#         i += 1
#         words = line.rstrip()
#         words = words.split("\t")
#         srna = words[-1]
#         start = words[2]
#         try:
#             start = int(start)
#         except ValueError:
#             continue
#         end = words[3]
#         end = int(end)
#         if end - start > 50:
#             strand = words[4]
#             new_feature = words[8]
#             feature = words[1]
#             overlap = words[7]
#             if new_feature == "FALSE":
#                 srna_type = "known"
#             else:
#                 srna_type = "novel"
#             srnaSeq = my_seq[start:end]
#             srnaSeqRev = srnaSeq.reverse_complement()
#             if strand == "-":
#                 srnaSeq = srnaSeqRev
#             if srna_type == "known":
#                 srnaPCFile = open("%s/positive_control/%s.fna" % (write_path, accession), "a")
#                 srnaPCFile.write(">%s[%s-%s,%s,%s,%s,%s]\n%s\n" % (srna, start, end, strand, srna_type, feature, overlap, srnaSeq))
#             else:
#                 srnaPredictedFile = open("%s/predicted/%s.fna" % (write_path, accession), "a")
#                 srnaPredictedFile.write(">%s[%s-%s,%s,%s,%s,%s]\n%s\n" % (srna, start, end, strand, srna_type, feature, overlap, srnaSeq))
def main():

    accession = rungetopts()

    try:
        inFile = open("/Users/thomasnicholson/phd/RNASeq/new_calls/%s_new_calls.txt" % accession, 'r')
    except IOError:
        print "/Users/thomasnicholson/phd/RNASeq/new_calls/%s_new_calls.txt not found" % accession
        sys.exit(2)
    try:
        fastaFile = list(SeqIO.parse("/Users/thomasnicholson/phd/RNASeq/sequences/%s.fna" % accession, "fasta"))
    except IOError:
        print "/Users/thomasnicholson/phd/RNASeq/sequences/%s.fna not found" % accession
        sys.exit(2)


    write_path = "/Users/thomasnicholson/phd/RNASeq/srna_seqs/version_1"

    print "Combining contigs"
    my_seq = srna.concatenateSequence(fastaFile)

    print "Writing sequences"
    srna.writeSequences(inFile,my_seq,accession,write_path)


if __name__ == "__main__":
    main()





#!/usr/bin/python

'''
file paths are hard coded
'''


import sys
from Bio import SeqIO
import getopt
import os
from BCBio import GFF
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import random
import comparativeSRNA as srna


help = '''
    {script_name} -c com_port [-o output_file] [--loglevel level]

    Reads the temperature data from a radio.  The temperature data is output in csv form.

    examples:
        Read table from radio attached to com4 and write the table to the file
        output.csv.

            {script_name} -c com4 -o output.csv

        Read table from radio attached to com3 and write the table to stdout. 
        You can use IO redirection to send the contents where every you want.

            # just print to the terminal 
            {script_name} -c com3

            # redirect to another file
            {script_name} -c com3 > somefile.csv

            # filter out temperatures that are -100
            {script_name} -c com3 | grep -v '^-100' 


    -c com_port
    --com_port comport
        Name of the COM port attached to the radio

    -o output_file
    --output output_file
        If specified write the table data to the given file.  If not specified
        the data will be written to stdout.

    --loglevel critical | error | warning | info | debug | notset
        Control the verbosity of the script by setting the log level.  Critical
        is the least verbose and notset is the most verbose.

        The default loglevel is {default_loglevel}.

        These values correspond directly to the python logging module levels. 
        (i.e. https://docs.python.org/3/howto/logging.html#logging-levels)


    -h 
    --help 
        print this message

'''

def usage():
    print help

def rungetopts():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "a:sqh", ["accession", "shuffle", "quiet", "help"])
    except getopt.GetoptError as err:
        # print help information and exit:
        print(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    accession = ""
    shuffled = False
    for o, a in opts:
            if o in ("-h", "--help"):
                usage()
                sys.exit()
            elif o in ("-a", "--accession"):
                accession = a
            elif o in ("-s", "--shuffle"):
                shuffled = True
            else:
                assert False, "unhandled option"
    if accession == "":
        print "-a <accession> missing. For more help use -h"
        sys.exit(2)
    return(accession, shuffled)



# def file_len(fname):
#     with open(fname) as f:
#         for i, l in enumerate(f):
#             pass
#     return i
#
# def intergenicSequence(accession, my_seq, shuffled):
#     start = 0
#     end = 0
#     random_seq = Seq("AG", generic_dna)
#     try:
#
#         in_handle = open("/Users/thomasnicholson/phd/RNASeq/sequences/%s.gff" % accession)
#         for rec in GFF.parse(in_handle):
#             for feature in rec.features:
#                 qualifiers = feature.qualifiers
#                 try:
#                     location = feature.location
#                     end = location.start
#                     intergeneicSeq = my_seq[start:end]
#                     if shuffled == True:
#                         shuffledSeq = ''.join(random.sample(str(intergeneicSeq), len(intergeneicSeq)))
#                         random_seq = random_seq + shuffledSeq
#                     else:
#                         random_seq = random_seq + intergeneicSeq
#                     start = location.end
#                     #print len(random_seq)
#                 except KeyError:
#                     pass
#
#         in_handle.close()
#
#     except IOError:
#         print "/Users/thomasnicholson/phd/RNASeq/sequences/%s.gff not found" % accession
#         sys.exit(2)
#     return random_seq
#
#
# def intergenicPositions(accession):
#     start = 0
#     end = 0
#     positions = [0]
#     try:
#
#         in_handle = open("/Users/thomasnicholson/phd/RNASeq/sequences/%s.gff" % accession)
#         for rec in GFF.parse(in_handle):
#             i = 0
#             for feature in rec.features:
#                 qualifiers = feature.qualifiers
#                 try:
#                     qualifiers['gene_biotype']
#                 except KeyError:
#                     continue
#                 try:
#                     i += 1
#                     location = feature.location
#                     end = location.start - 49
#                     if end < start:
#                         continue
#                     tmpPos = range(start,end)
#                     positions = positions + tmpPos
#                     start = location.end + 50
#                 except KeyError:
#                     pass
#
#         in_handle.close()
#
#     except IOError:
#         print "/Users/thomasnicholson/phd/RNASeq/sequences/%s.gff not found" % accession
#         sys.exit(2)
#     return positions
#
#
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
# def selectRandomLocation(inFile, positions,fileLength, random_seq, accession):
#
#     randomFile = open("/Users/thomasnicholson/phd/RNASeq/new_calls/random/python_version_1/%s_random_no_shuffle_new_calls.txt" % accession, "w")
#     randomFile.write("start\tend\tstrand\tsequence\n")
#
#     shuffledIndexes = random.sample(positions, fileLength)
#     seqLength = len(random_seq)
#     seqIndexes = random.sample(range(0,seqLength), fileLength)
#
#     srnaLengths = []
#     srnaStrands = []
#     srnaIDs = []
#     i = 0
#     for line in inFile:
#         i += 1
#         words = line.rstrip()
#         words = words.split("\t")
#         start = words[2]
#         try:
#             start = int(start)
#         except ValueError:
#             continue
#         end = words[3]
#         end = int(end)
#         srna = words[-1]
#         srna_length = end - start
#         srnaLengths.append(srna_length)
#         strand = words[4]
#         srnaStrands.append(strand)
#         srnaIDs.append(srna)
#     for i in range(0,len(shuffledIndexes)):
#         index = shuffledIndexes[i]
#         length = srnaLengths[i]
#         strand = srnaStrands[i]
#         seqIndex = seqIndexes[i]
#         srna = srnaIDs[i]
#         if strand == "+":
#             start = index
#             end = start + length
#             seqStart = seqIndex
#             seqEnd  = seqStart + length
#         else:
#             end = index
#             start = end - length
#             seqEnd = seqIndex
#             seqStart  = seqEnd - length
#         if start < 1:
#             continue
#         if end < 1:
#             continue
#         if length < 50:
#             continue
#         if length > 500:
#             continue
#         if seqStart < 1:
#             continue
#         if seqEnd < 1:
#             continue
#         if seqEnd > seqLength:
#             continue
#         if seqStart > seqLength:
#             continue
#         sequence  = random_seq[seqStart:seqEnd]
#         randomFile = open("/Users/thomasnicholson/phd/RNASeq/new_calls/random/python_version_1/%s_random_no_shuffle_new_calls.txt" % accession, "a")
#         randomFile.write("%s\t%s\t%s\t%s\n" % (start, end, strand, sequence))
#
#         srna_type = "random"
#
#         write_path = "/Users/thomasnicholson/phd/RNASeq/srna_seqs/version_1/negative_control_no_shuffle"
#         srnaFile = open("%s/%s.fna" % (write_path, accession), "a")
#         srnaFile.write(">%s[%s-%s,%s,%s]\n%s\n" % (srna, seqStart, seqEnd, strand, srna_type, sequence))

def main():

    accession, shuffled = rungetopts()
    print "Reading files"
    try:
        inFile = open("/Users/thomasnicholson/phd/RNASeq/new_calls/%s_new_calls.txt" % accession, 'r')

        fileLength = file_len("/Users/thomasnicholson/phd/RNASeq/new_calls/%s_new_calls.txt" % accession)
    except IOError:
        print "/Users/thomasnicholson/phd/RNASeq/new_calls/%s_new_calls.txt not found" % accession
        sys.exit(2)
    try:
        fastaFile = list(SeqIO.parse("/Users/thomasnicholson/phd/RNASeq/sequences/%s.fna" % accession, "fasta"))
    except IOError:
        print "/Users/thomasnicholson/phd/RNASeq/sequences/%s.fna not found" % accession
        sys.exit(2)

    print "Combining contigs"
    my_seq = srna.concatenateSequence(fastaFile)



    print "Getting intergenic sequence"
    random_seq = srna.intergenicSequence(accession, my_seq, shuffled)


    print "Getting intergenic positions"
    positions = srna.intergenicPositions(accession)

    print "Selecting random sRNAs"
    srna.selectRandomLocation(inFile, positions,fileLength, random_seq, accession)





if __name__ == "__main__":
    main()





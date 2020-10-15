#!/usr/bin/python
import sys
from Bio import SeqIO
import getopt
import os
from BCBio import GFF
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import random
import pandas as pd
import comparativesrna as srna

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
        opts, args = getopt.getopt(sys.argv[1:], "a:rqh", ["accession", "random", "quiet", "help"])
    except getopt.GetoptError as err:
        # print help information and exit:
        print(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    accession = ""
    random = False
    for o, a in opts:
            if o in ("-h", "--help"):
                usage()
                sys.exit()
            elif o in ("-a", "--accession"):
                accession = a
            elif o in ("-r", "--random"):
                random = True
            else:
                assert False, "unhandled option"
    if accession == "":
        print "-a <accession> missing. For more help use -h"
        sys.exit(2)
    return(accession, random)

# def getreaddepths(accession):
#     try:
#         df = None
#         for filename in os.listdir("/Users/thomasnicholson/phd/RNASeq/plot_files/%s/" % accession):
#             filesize = os.path.getsize(
#                 "/Users/thomasnicholson/phd/RNASeq/plot_files/%s/%s" % (accession, filename))
#             if filesize == 0:
#                 print("No data in %s" % filename)
#                 continue
#             plotFile = pd.read_csv(
#                 os.path.join("/Users/thomasnicholson/phd/RNASeq/plot_files/%s/" % accession, filename),
#                 sep='\t', header=None)
#             print(filename)
#             plotFile['selected'] = plotFile.iloc[:].max(axis=1)
#             tmpDf = plotFile.iloc[:, 2]
#             if df is not None:
#                 df = pd.concat([df.reset_index(drop=True), tmpDf], axis=1)
#             else:
#                 df = tmpDf
#         dfOut = df
#         dfOut['mean'] = df.iloc[:].mean(axis=1)
#         dfOut['median'] = df.median(axis=1)
#         dfOut['max'] = df.max(axis=1)
#         return dfOut
#
#     except IOError:
#         print "Cannot open a file in /Users/thomasnicholson/phd/RNASeq/plot_files/%s/" % accession
#
# def sRNA_read_depths(inFile, read_depths_df,accession, random):
#     if random == False:
#         outFile  = open("/Users/thomasnicholson/phd/RNASeq/srna_seqs/version_1/read_depths/%s_read_depths.txt" % accession, 'w')
#         outFile.write("ID\tstart\tend\tgroup\tfeature\tmean_mean\tmean_median\tmean_max\tmedian_mean\tmedian_median\tmedian_max\tmax_mean\tmax_median\tmax_max\n")
#         outFile.close()
#         outFile  = open("/Users/thomasnicholson/phd/RNASeq/srna_seqs/version_1/read_depths/%s_read_depths.txt" % accession, 'a')
#     else:
#         outFile  = open("/Users/thomasnicholson/phd/RNASeq/srna_seqs/version_1/read_depths_negative_control/%s_read_depths.txt" % accession, 'w')
#         outFile.write("ID\tstart\tend\tgroup\tfeature\tmean_mean\tmean_median\tmean_max\tmedian_mean\tmedian_median\tmedian_max\tmax_mean\tmax_median\tmax_max\n")
#         outFile.close()
#         outFile  = open("/Users/thomasnicholson/phd/RNASeq/srna_seqs/version_1/read_depths_negative_control/%s_read_depths.txt" % accession, 'a')
#
#     if random == False:
#         for line in inFile:
#             words = line.rstrip()
#             words = words.split("\t")
#             srna = words[-1]
#             start = words[2]
#             try:
#                 start = int(start)
#             except ValueError:
#                 continue
#             end = words[3]
#             end = int(end)
#             new_feature = words[8]
#             feature = words[1]
#             if new_feature == "FALSE":
#                 srna_type = "known"
#             else:
#                 srna_type = "novel"
#
#             subsetDF = read_depths_df[start:end]
#
#             mean_mean = subsetDF['mean'].mean()
#             mean_median = subsetDF['mean'].median()
#             mean_max = subsetDF['mean'].max()
#             median_mean = subsetDF['median'].mean()
#             median_median = subsetDF['median'].median()
#             median_max = subsetDF['median'].mean()
#             max_mean = subsetDF['max'].mean()
#             max_median = subsetDF['max'].median()
#             max_max = subsetDF['max'].max()
#
#             outFile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (srna, start, end, srna_type, feature, mean_mean, mean_median, mean_max, median_mean, median_median, median_max, max_mean, max_median, max_max))
#     else:
#         i = 0
#         for line in inFile:
#             i += 1
#             words = line.rstrip()
#             words = words.split("\t")
#             srna = "%s_%s" % (accession, i)
#             start = words[0]
#             try:
#                 start = int(start)
#             except ValueError:
#                 continue
#             end = words[1]
#             end = int(end)
#             feature = "intergenic"
#             srna_type = "negative_control"
#
#             subsetDF = read_depths_df[start:end]
#
#             mean_mean = subsetDF['mean'].mean()
#             mean_median = subsetDF['mean'].median()
#             mean_max = subsetDF['mean'].max()
#             median_mean = subsetDF['median'].mean()
#             median_median = subsetDF['median'].median()
#             median_max = subsetDF['median'].mean()
#             max_mean = subsetDF['max'].mean()
#             max_median = subsetDF['max'].median()
#             max_max = subsetDF['max'].max()
#
#             outFile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
#             srna, start, end, srna_type, feature, mean_mean, mean_median, mean_max, median_mean, median_median,
#             median_max, max_mean, max_median, max_max))



def main():

    accession, random = rungetopts()
    print "Reading files"
    try:
        if random == False:
            inFile = open("/Users/thomasnicholson/phd/RNASeq/new_calls/%s_new_calls.txt" % accession, 'r')
        else:
            inFile = open("/Users/thomasnicholson/phd/RNASeq/new_calls/random/python_version_1/%s_random_new_calls.txt" % accession, 'r')
    except IOError:
        print "/Users/thomasnicholson/phd/RNASeq/new_calls/%s_new_calls.txt not found" % accession
        sys.exit(2)

    print "getting read depths"
    read_depths_df = srna.getreaddepths(accession)
    # read_depths_df.to_csv("/Users/thomasnicholson/phd/RNASeq/srna_seqs/version_1/read_depths/%s_read_depths.csv" % accession)

    print "writing file"

    srna.sRNA_read_depths(inFile, read_depths_df, accession, random)


if __name__ == "__main__":
    main()



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
        opts, args = getopt.getopt(sys.argv[1:], "a:qh", ["accession", "quiet", "help"])
    except getopt.GetoptError as err:
        # print help information and exit:
        print(err) # will print something like "option -a not recognized"
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

def getreaddepths(accession):
    try:
        read_depths_list = {}
        for filename in os.listdir("/Users/thomasnicholson/phd/RNASeq/plot_files/%s/" % accession):
                plotFile  = open(os.path.join("/Users/thomasnicholson/phd/RNASeq/plot_files/%s/" % accession, filename), 'r')
                print filename

                read_depths = []
                for line in plotFile:
                    words = line.rstrip()
                    words = words.split()
                    pos = words[0]
                    neg = words[1]
                    if pos > neg:
                        read_depths.append(pos)
                    else:
                        read_depths.append(neg)
                read_depths_list[filename] = read_depths

        df = pd.DataFrame(data = read_depths_list)
        df['mean'] = df.mean(axis=1)
        df['median'] = df.median(axis=1)
        df['max'] = df.max(axis=1)
        return df
    except IOError:
        print "Cannot open a file in /Users/thomasnicholson/phd/RNASeq/plot_files/%s/" % accession

def sRNA_read_depths(inFile, read_depths_df,accession):
    outFile  = open("/Users/thomasnicholson/phd/RNASeq/srna_seqs/version_1/read_depths/%s_read_depths.txt" % accession, 'w')
    outFile.write("ID\tstart\tend\tgroup\tfeature\tmean_mean\tmean_median\tmean_max\tmedian_mean\tmedian_median\tmedian_max\tmax_mean\tmax_median\tmax_max\n")
    outFile.close()
    outFile  = open("/Users/thomasnicholson/phd/RNASeq/srna_seqs/version_1/read_depths/%s_read_depths.txt" % accession, 'a')

    for line in inFile:
        words = line.rstrip()
        words = words.split("\t")
        srna = words[-1]
        start = words[2]
        try:
            start = int(start)
        except ValueError:
            continue
        end = words[3]
        end = int(end)
        new_feature = words[8]
        feature = words[1]
        if new_feature == "FALSE":
            srna_type = "known"
        else:
            srna_type = "novel"

        subsetDF = read_depths_df[start:end]

        mean_mean = subsetDF['mean'].mean()
        mean_median = subsetDF['mean'].median()
        mean_max = subsetDF['mean'].max()
        median_mean = subsetDF['median'].mean()
        median_median = subsetDF['median'].median()
        median_max = subsetDF['median'].mean()
        max_mean = subsetDF['max'].mean()
        max_median = subsetDF['max'].median()
        max_max = subsetDF['max'].max()

        outFile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (srna, start, end, srna_type, feature, mean_mean, mean_median, mean_max, median_mean, median_median, median_max, max_mean, max_median, max_max))




def main():

    accession = rungetopts()
    print "Reading files"
    try:
        inFile = open("/Users/thomasnicholson/phd/RNASeq/new_calls/%s_new_calls.txt" % accession, 'r')

    except IOError:
        print "/Users/thomasnicholson/phd/RNASeq/new_calls/%s_new_calls.txt not found" % accession
        sys.exit(2)

    read_depths_df = getreaddepths(accession)

    sRNA_read_depths(inFile, read_depths_df, accession)


if __name__ == "__main__":
    main()



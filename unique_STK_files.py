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
        opts, args = getopt.getopt(sys.argv[1:], "f:qh", ["file", "quiet", "help"])
    except getopt.GetoptError as err:
        # print help information and exit:
        print(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    file = ""
    for o, a in opts:
            if o in ("-h", "--help"):
                usage()
                sys.exit()
            elif o in ("-f", "--file"):
                file = a
            else:
                assert False, "unhandled option"
    if file == "":
        print "-f <file> missing. For more help use -h"
        sys.exit(2)
    return(file)

def unique_stk(inFile):
    new = False
    writeFile = True
    list = []
    for line in inFile:

        words = line.rstrip()
        if words == "# STOCKHOLM 1.0":
            new = True
            writeFile = True
            try:
                outFile.close()
            except UnboundLocalError:
                pass
            continue
        if new == True and "#=GS" in line:
            names = words.split(" ")
            id = names[-1]
            id = id.split("[")
            id = id[0]
            subsetID = names[1]
            if "|" in subsetID:
                subsetID = subsetID.split("|")
                subsetID = subsetID[1]
            if subsetID in list:
                writeFile = False
                continue
            list.append(subsetID)
            print id
            print subsetID
            if writeFile == False:
                continue
            outFile = open("/Users/thomasnicholson/phd/RNASeq/srna_seqs/version_1/stockholm/%s.stk" % id, 'a')
            outFile.write("# STOCKHOLM 1.0\n\n%s\n" % words)
            new = False
        else:
            try:
                outFile.write("%s\n" % words)
            except UnboundLocalError:
                pass
            except ValueError:
                pass


def main():

    filename = rungetopts()
    print "Reading files"
    try:
        inFile = open("/Users/thomasnicholson/phd/RNASeq/srna_seqs/version_1/%s" % filename, 'r')

    except IOError:
        print "/Users/thomasnicholson/phd/RNASeq/srna_seqs/version_1/%s not found" % filename
        sys.exit(2)

    unique_stk(inFile)


if __name__ == "__main__":
    main()
#!/usr/bin/python
import sys
from Bio import SeqIO
from Bio.Alphabet import generic_dna
import getopt
import os
# args = '-c -f -p -q -h'.split()
# optlist, args = getopt.getopt(args, 'c:f:p:qh')


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

def main():
    input_name = ""
    fasta_name = ""
    file_path = "."
    try:
        opts, args = getopt.getopt(sys.argv[1:], "a:c:f:p:qh", ["accession", "calls", "fasta", "filepath", "quiet", "help"])
    except getopt.GetoptError as err:
        # print help information and exit:
        print(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    for o, a in opts:
            if o in ("-h", "--help"):
                usage()
                sys.exit()
            elif o in ("-c", "--calls"):
                input_name = a
            elif o in ("-a", "--accession"):
                accession = a
            elif o in ("-f", "--fasta"):
                fasta_name = a
            elif o in ("-p", "--filepath"):
                file_path = a
            else:
                assert False, "unhandled option"

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
        
    my_seq = fastaFile[0].seq
    write_path = "/Users/thomasnicholson/phd/RNASeq/srna_seqs/version_1/%s/" % accession
    if os.path.isdir(write_path) == False:
        try:
            os.mkdir(write_path)
        except OSError:
            print ("Creation of the directory %s failed" % accession)
            sys.exit(2)
    directory = os.listdir(write_path)
    if len(directory) != 0:
        print "Examples of files in %s" % write_path
        print directory[0:4]
        query_user = raw_input("%s is not an empty directory. Continue anyway y/n (this may write over existing files): " % write_path)
        if query_user == "y":
            print "Using %s as directory" % write_path
        else:
            print "Exiting script"
            sys.exit(2)


    i = 0
    for seq in fastaFile:
        if i == 0:
            i += 1
            continue
        i += 1
        my_seq = my_seq + seq.seq
    i = 0
    for line in inFile:
        i += 1
        words = line.rstrip()
        words = words.split("\t")
        id = words[0]
        srna = words[-1]
        start = words[2]
        try:
            start = int(start)
        except ValueError:
            continue
        end = words[3]
        end = int(end)
        if end - start > 50:
            strand = words[4]
            new_feature = words[8]
            feature = words[1]
            overlap = words[7]
            if new_feature == "FALSE":
                srna_type = "known"
            else:
                srna_type = "novel"
            srnaSeq = my_seq[start:end]
            srnaSeqRev = srnaSeq.reverse_complement()
            if strand == "-":
                srnaSeq = srnaSeqRev
            srnaFile = open("%s/%s.fna" % (write_path, srna), "a")
            srnaFile.write(">%s[%s-%s,%s,%s,%s,%s]\n%s\n" %(srna,start, end, strand, srna_type, feature, overlap, srnaSeq))


if __name__ == "__main__":
    main()





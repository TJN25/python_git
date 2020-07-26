#!/usr/bin/python

from Bio import SeqIO

from BCBio import GFF
import sys

try:
    sys.argv[1]
except IndexError:
    print "Error: A genbank file needs to be included in the argument parameters"
    exit()

names = sys.argv[1]
print(names)
names = names.split(".")[0]
print( "%s.gff" % names)

try:
    output_name = sys.argv[2]
except IndexError:
    names = sys.argv[1]
    names = names.split(".")[0]
    output_name = "%s.gff3" % names



in_file = sys.argv[1]
out_file = output_name
in_handle = open(in_file)
out_handle = open(out_file, "w")
print "File read in: %s" % sys.argv[1]

GFF.write(SeqIO.parse(in_handle, "embl"), out_handle)

in_handle.close()
out_handle.close()




print "Writing output to: %s" % output_name

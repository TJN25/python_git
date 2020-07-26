#!/usr/bin/python

from Bio import SeqIO

from BCBio import GFF



in_file = "/Users/thomasnicholson/Downloads/KJ183193.1.txt"
out_file = "/Users/thomasnicholson/Downloads/KJ183193.1.gff"
in_handle = open(in_file)
out_handle = open(out_file, "w")

GFF.write(SeqIO.parse(in_handle, "embl"), out_handle)

in_handle.close()
out_handle.close()
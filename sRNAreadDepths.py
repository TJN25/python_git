#!/usr/bin/python
import sys
import getopt

import pandas as pd
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
        opts, args = getopt.getopt(sys.argv[1:], "n:r:o:qh", ["nhmmer", "readdepths", "output", "quiet", "help"])
    except getopt.GetoptError as err:
        # print help information and exit:
        print(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    nhmmer = ""
    readdepths = ""
    output = ""

    for o, a in opts:
            if o in ("-h", "--help"):
                usage()
                sys.exit()
            elif o in ("-n", "--nhmmer"):
                nhmmer = a
            elif o in ("-r", "--readdepths"):
                readdepths = a
            elif o in ("-o", "--output"):
                output = a
            else:
                assert False, "unhandled option"
    if nhmmer == "":
        print "-n <nhmmer> missing. For more help use -h"
        sys.exit(2)
    if readdepths == "":
        print "-r <readdepths> missing. For more help use -h"
        sys.exit(2)
    if output == "":
        print "-o <output> missing. For more help use -h"
        sys.exit(2)
    return(nhmmer, readdepths, output)


# def openNHMMER(nhmmername):
#     nhmmerDF = pd.read_csv(nhmmername, delim_whitespace=True, header=None, comment='#')
#     nhmmerDF.columns = ["target_name", "accession", "query_name", "accession_2", "hmmfrom", "hmmto", "alifrom", "alito", "envfrom", "envto", "sq_len", "strand", "E_value", "score", "bias", "description_of_target"]
#     nhmmerDF[["ID", "descriptors"]] = nhmmerDF.target_name.str.split("[", expand = True)
#     nhmmerDF[["ID_2", "descriptors_2"]] = nhmmerDF.query_name.str.split("[", expand = True)
#     d = nhmmerDF.groupby('ID')['ID_2'].apply(list).to_dict()
#     return(d)
#
# def openReadDepths(readdepthsname, d):
#     readdepthsDF = pd.read_csv(readdepthsname, sep = "\t", comment='#')
#     readdepthsDF = readdepthsDF[readdepthsDF['ID'] != "ID"]
#
#     ##when being done in jupyter the columns were all read in as string and the lines below were necessary...
#     ##it seems to work fine now
#
#     # print(readdepthsDF.dtypes)
#     # readdepthsDF[["mean_value", "mean_decimal"]] = readdepthsDF.max_mean.str.split(".", expand = True)
#     # print(1)
#     # readdepthsDF[["median_value", "median_decimal"]] = readdepthsDF.max_median.str.split(".", expand = True)
#     # readdepthsDF[["max_value", "max_decimal"]] = readdepthsDF.max_max.str.split(".", expand = True)
#     # readdepthsDF[['mean_value', 'median_value', 'max_value']] = readdepthsDF.loc[:,['mean_value', 'median_value', 'max_value']].apply(pd.to_numeric)
#
#
#     readdepthsDF["mean_value"] = readdepthsDF['max_mean']
#     readdepthsDF["median_value"] = readdepthsDF['max_median']
#     readdepthsDF["max_value"] = readdepthsDF['max_max']
#     readdepthsDF[['mean_value', 'median_value', 'max_value']] = readdepthsDF.loc[:,['mean_value', 'median_value', 'max_value']].apply(pd.to_numeric)
#
#
#     idList = list(d.keys())
#     readdepthsKept = readdepthsDF[readdepthsDF['ID'].isin(idList)]
#     return(readdepthsKept)
#
#
# def writeReadDepths(outname, readDepths, d):
#     seen = []
#     d2 = {}
#     i = 0
#
#     outFile = open(outname, "w")
#     outFile.write(
#         "ID\tmean_mean\tmean_median\tmean_max\tmedian_mean\tmedian_median\tmedian_max\tmax_mean\tmax_median\tmax_max\tID_2\n")
#     outFile.close()
#     outFile = open(outname, "a")
#     values = []
#     for key in d:
#         #     print(i)
#         #     i += 1
#         #     if i > 100:
#         #         break
#         #     if key in seen:
#         #         continue
#         #     print(key)
#         #     print(seen)
#         values = d[key]
#         seen.append(values)
#         df = readDepths[readDepths['ID'].isin(values)]
#         #     print(df['mean_value'].dtypes)
#
#         #     print(df['mean_value'].dtypes)
#         mean_mean = df['mean_value'].mean()
#         mean_median = df['mean_value'].median()
#         mean_max = df['mean_value'].max()
#         median_mean = df['median_value'].mean()
#         median_median = df['median_value'].median()
#         median_max = df['median_value'].max()
#         max_mean = df['max_value'].mean()
#         max_median = df['max_value'].median()
#         max_max = df['max_value'].max()
#         #     print(key)
#         #     print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (key,mean_mean,mean_median,mean_max,median_mean,median_median,median_max,max_mean,max_median,max_max,values))
#         outFile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
#         key, mean_mean, mean_median, mean_max, median_mean, median_median, median_max, max_mean, max_median, max_max,
#         values))
#     outFile.close()

def main():

    nhmmer, readdepths, output = rungetopts()

    print("Opening nhmmer results: %s" % nhmmer)
    d = srna.openNHMMER(nhmmername = nhmmer)

    print("Opening read depths results: %s" % readdepths)
    readDepthsKept = srna.openReadDepths(readdepthsname=readdepths, d=d)

    print("Writing file: %s" % output)
    srna.writeReadDepths(outname = output, readDepths = readDepthsKept, d = d)


if __name__ == "__main__":
    main()


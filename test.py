#!/usr/bin/python

import comparativeSRNA as srna


help = '''
    foobar.py v 0.1 (August 2020) is a script for {}.
    Usage:
         foobar.py [options] <input> <output>
    
    Options:
        -h	Display this help
        -q  Supress messages
    Input
        -i	<input> the input
        -o	<output> the output

    
'''

def usage():
    print help

def rungetopts():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "i:o:qh", ["input", "output", "quiet", "help"])
    except getopt.GetoptError as err:
        print(err)
        usage()
        sys.exit(2)
    input = ""
    output = ""
    for o, a in opts:
            if o in ("-h", "--help"):
                usage()
                sys.exit()
            elif o in ("-i", "--input"):
                input = a
            elif o in ("-o", "--output"):
                output = a
            else:
                assert False, "unhandled option"
    if output == "":
        print "-o <output> missing. For more help use -h"
        sys.exit(2)
    if input == "":
        print "-i <input> missing. For more help use -h"
        sys.exit(2)
    return(input, output)


srna.helloworld()
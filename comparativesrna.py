import os
import sys
import random

import pandas as pd
from BCBio import GFF
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

def package_test():
	print("comparativesrna.py loaded")

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i


def intergenicSequence(accession, my_seq, shuffled):
    start = 0
    end = 0
    random_seq = Seq("AG", generic_dna)
    try:

        in_handle = open("/Users/thomasnicholson/phd/RNASeq/sequences/%s.gff" % accession)
        for rec in GFF.parse(in_handle):
            for feature in rec.features:
                qualifiers = feature.qualifiers
                try:
                    location = feature.location
                    end = location.start
                    intergeneicSeq = my_seq[start:end]
                    if shuffled == True:
                        shuffledSeq = ''.join(random.sample(str(intergeneicSeq), len(intergeneicSeq)))
                        random_seq = random_seq + shuffledSeq
                    else:
                        random_seq = random_seq + intergeneicSeq
                    start = location.end
                    #print(len(random_seq))
                except KeyError:
                    pass

        in_handle.close()

    except IOError:
        print("/Users/thomasnicholson/phd/RNASeq/sequences/%s.gff not found" % accession)
        sys.exit(2)
    return random_seq


def intergenicPositions(accession):
    start = 0
    end = 0
    positions = [0]
    try:

        in_handle = open("/Users/thomasnicholson/phd/RNASeq/sequences/%s.gff" % accession)
        for rec in GFF.parse(in_handle):
            i = 0
            for feature in rec.features:
                qualifiers = feature.qualifiers
                try:
                    qualifiers['gene_biotype']
                except KeyError:
                    continue
                try:
                    i += 1
                    location = feature.location
                    end = location.start - 49
                    if end < start:
                        continue
                    tmpPos = range(start,end)
                    positions = positions + tmpPos
                    start = location.end + 50
                except KeyError:
                    pass

        in_handle.close()

    except IOError:
        print("/Users/thomasnicholson/phd/RNASeq/sequences/%s.gff not found" % accession)
        sys.exit(2)
    return positions


def makeoutputdirectory(write_path):
    if os.path.isdir(write_path) == False:
        try:
            os.mkdir(write_path)
        except OSError:
            print("Creation of the directory %s failed" % write_path)
            sys.exit(2)
    directory = os.listdir(write_path)
    if len(directory) != 0:
        print("Examples of files in %s" % write_path)
        print(directory[0:4])
        query_user = input("%s is not an empty directory. Continue anyway y/n (this may write over existing files): " % write_path)
        if query_user == "y":
            print("Using %s as directory" % write_path)
        else:
            print("Exiting script")
            sys.exit(2)


def concatenateSequence(fastaFile):
    my_seq = fastaFile[0].seq
    i = 0
    for seq in fastaFile:
        if i == 0:
            i += 1
            continue
        i += 1
        my_seq = my_seq + seq.seq
    return my_seq


def selectRandomLocation(inFile, positions,fileLength, random_seq, accession):

    randomFile = open("/Users/thomasnicholson/phd/RNASeq/new_calls/random/python_version_1/%s_random_no_shuffle_new_calls.txt" % accession, "w")
    randomFile.write("start\tend\tstrand\tsequence\n")

    shuffledIndexes = random.sample(positions, fileLength)
    seqLength = len(random_seq)
    seqIndexes = random.sample(range(0,seqLength), fileLength)

    srnaLengths = []
    srnaStrands = []
    srnaIDs = []
    i = 0
    for line in inFile:
        i += 1
        words = line.rstrip()
        words = words.split("\t")
        start = words[2]
        try:
            start = int(start)
        except ValueError:
            continue
        end = words[3]
        end = int(end)
        srna = words[-1]
        srna_length = end - start
        srnaLengths.append(srna_length)
        strand = words[4]
        srnaStrands.append(strand)
        srnaIDs.append(srna)
    for i in range(0,len(shuffledIndexes)):
        index = shuffledIndexes[i]
        length = srnaLengths[i]
        strand = srnaStrands[i]
        seqIndex = seqIndexes[i]
        srna = srnaIDs[i]
        if strand == "+":
            start = index
            end = start + length
            seqStart = seqIndex
            seqEnd  = seqStart + length
        else:
            end = index
            start = end - length
            seqEnd = seqIndex
            seqStart  = seqEnd - length
        if start < 1:
            continue
        if end < 1:
            continue
        if length < 50:
            continue
        if length > 500:
            continue
        if seqStart < 1:
            continue
        if seqEnd < 1:
            continue
        if seqEnd > seqLength:
            continue
        if seqStart > seqLength:
            continue
        sequence  = random_seq[seqStart:seqEnd]
        randomFile = open("/Users/thomasnicholson/phd/RNASeq/new_calls/random/python_version_1/%s_random_no_shuffle_new_calls.txt" % accession, "a")
        randomFile.write("%s\t%s\t%s\t%s\n" % (start, end, strand, sequence))

        srna_type = "random"

        write_path = "/Users/thomasnicholson/phd/RNASeq/srna_seqs/version_1/negative_control_no_shuffle"
        srnaFile = open("%s/%s.fna" % (write_path, accession), "a")
        srnaFile.write(">%s[%s-%s,%s,%s]\n%s\n" % (srna, seqStart, seqEnd, strand, srna_type, sequence))


def getreaddepths(accession):
    try:
        df = None
        for filename in os.listdir("/Users/thomasnicholson/phd/RNASeq/plot_files/%s/" % accession):
            filesize = os.path.getsize(
                "/Users/thomasnicholson/phd/RNASeq/plot_files/%s/%s" % (accession, filename))
            if filesize == 0:
                print("No data in %s" % filename)
                continue
            plotFile = pd.read_csv(
                os.path.join("/Users/thomasnicholson/phd/RNASeq/plot_files/%s/" % accession, filename),
                sep='\t', header=None)
            print(filename)
            plotFile['selected'] = plotFile.iloc[:].max(axis=1)
            tmpDf = plotFile.iloc[:, 2]
            if df is not None:
                df = pd.concat([df.reset_index(drop=True), tmpDf], axis=1)
            else:
                df = tmpDf
        dfOut = df
        dfOut['mean'] = df.iloc[:].mean(axis=1)
        dfOut['median'] = df.median(axis=1)
        dfOut['max'] = df.max(axis=1)
        return dfOut

    except IOError:
        print("Cannot open a file in /Users/thomasnicholson/phd/RNASeq/plot_files/%s/" % accession)


def sRNA_read_depths(inFile, read_depths_df,accession, random):
    if random == False:
        outFile  = open("/Users/thomasnicholson/phd/RNASeq/srna_seqs/version_1/read_depths/%s_read_depths.txt" % accession, 'w')
        outFile.write("ID\tstart\tend\tgroup\tfeature\tmean_mean\tmean_median\tmean_max\tmedian_mean\tmedian_median\tmedian_max\tmax_mean\tmax_median\tmax_max\n")
        outFile.close()
        outFile  = open("/Users/thomasnicholson/phd/RNASeq/srna_seqs/version_1/read_depths/%s_read_depths.txt" % accession, 'a')
    else:
        outFile  = open("/Users/thomasnicholson/phd/RNASeq/srna_seqs/version_1/read_depths_negative_control/%s_read_depths.txt" % accession, 'w')
        outFile.write("ID\tstart\tend\tgroup\tfeature\tmean_mean\tmean_median\tmean_max\tmedian_mean\tmedian_median\tmedian_max\tmax_mean\tmax_median\tmax_max\n")
        outFile.close()
        outFile  = open("/Users/thomasnicholson/phd/RNASeq/srna_seqs/version_1/read_depths_negative_control/%s_read_depths.txt" % accession, 'a')

    if random == False:
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
    else:
        i = 0
        for line in inFile:
            i += 1
            words = line.rstrip()
            words = words.split("\t")
            srna = "%s_%s" % (accession, i)
            start = words[0]
            try:
                start = int(start)
            except ValueError:
                continue
            end = words[1]
            end = int(end)
            feature = "intergenic"
            srna_type = "negative_control"

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

            outFile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
            srna, start, end, srna_type, feature, mean_mean, mean_median, mean_max, median_mean, median_median,
            median_max, max_mean, max_median, max_max))


def single_fasta(fastaFile, folder):
    for seq in fastaFile:
        id = seq.id
        outname = id.split("[")
        outname = outname[0]
        my_seq = seq.seq
        outFile = open("/Users/thomasnicholson/phd/RNASeq/srna_seqs/version_1/%s/%s.fna" % (folder, outname), "w")
        outFile.write(">%s\n%s\n" % (id, my_seq))


def openNHMMER(nhmmername):
    nhmmerDF = pd.read_csv(nhmmername, delim_whitespace=True, header=None, comment='#')
    nhmmerDF.columns = ["target_name", "accession", "query_name", "accession_2", "hmmfrom", "hmmto", "alifrom", "alito", "envfrom", "envto", "sq_len", "strand", "E_value", "score", "bias", "description_of_target"]
    nhmmerDF[["ID", "descriptors"]] = nhmmerDF.target_name.str.split("[", expand = True)
    nhmmerDF[["ID_2", "descriptors_2"]] = nhmmerDF.query_name.str.split("[", expand = True)
    d = nhmmerDF.groupby('ID')['ID_2'].apply(list).to_dict()
    return(d)


def openReadDepths(readdepthsname, d):
    readdepthsDF = pd.read_csv(readdepthsname, sep = "\t", comment='#')
    readdepthsDF = readdepthsDF[readdepthsDF['ID'] != "ID"]

    ##when being done in jupyter the columns were all read in as string and the lines below were necessary...
    ##it seems to work fine now

    # print(readdepthsDF.dtypes)
    # readdepthsDF[["mean_value", "mean_decimal"]] = readdepthsDF.max_mean.str.split(".", expand = True)
    # print(1)
    # readdepthsDF[["median_value", "median_decimal"]] = readdepthsDF.max_median.str.split(".", expand = True)
    # readdepthsDF[["max_value", "max_decimal"]] = readdepthsDF.max_max.str.split(".", expand = True)
    # readdepthsDF[['mean_value', 'median_value', 'max_value']] = readdepthsDF.loc[:,['mean_value', 'median_value', 'max_value']].apply(pd.to_numeric)


    readdepthsDF["mean_value"] = readdepthsDF['max_mean']
    readdepthsDF["median_value"] = readdepthsDF['max_median']
    readdepthsDF["max_value"] = readdepthsDF['max_max']
    readdepthsDF[['mean_value', 'median_value', 'max_value']] = readdepthsDF.loc[:,['mean_value', 'median_value', 'max_value']].apply(pd.to_numeric)


    idList = list(d.keys())
    readdepthsKept = readdepthsDF[readdepthsDF['ID'].isin(idList)]
    return(readdepthsKept)


def writeReadDepths(outname, readDepths, d):
    seen = []
    d2 = {}
    i = 0

    outFile = open(outname, "w")
    outFile.write(
        "ID\tmean_mean\tmean_median\tmean_max\tmedian_mean\tmedian_median\tmedian_max\tmax_mean\tmax_median\tmax_max\tID_2\n")
    outFile.close()
    outFile = open(outname, "a")
    values = []
    for key in d:
        #     print(i)
        #     i += 1
        #     if i > 100:
        #         break
        #     if key in seen:
        #         continue
        #     print(key)
        #     print(seen)
        values = d[key]
        seen.append(values)
        df = readDepths[readDepths['ID'].isin(values)]
        #     print(df['mean_value'].dtypes)

        #     print(df['mean_value'].dtypes)
        mean_mean = df['mean_value'].mean()
        mean_median = df['mean_value'].median()
        mean_max = df['mean_value'].max()
        median_mean = df['median_value'].mean()
        median_median = df['median_value'].median()
        median_max = df['median_value'].max()
        max_mean = df['max_value'].mean()
        max_median = df['max_value'].median()
        max_max = df['max_value'].max()
        #     print(key)
        #     print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (key,mean_mean,mean_median,mean_max,median_mean,median_median,median_max,max_mean,max_median,max_max,values))
        outFile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
        key, mean_mean, mean_median, mean_max, median_mean, median_median, median_max, max_mean, max_median, max_max,
        values))
    outFile.close()


def writeSequences(inFile,my_seq,accession,write_path):
    i = 0
    for line in inFile:
        i += 1
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
            if srna_type == "known":
                srnaPCFile = open("%s/positive_control/%s.fna" % (write_path, accession), "a")
                srnaPCFile.write(">%s[%s-%s,%s,%s,%s,%s]\n%s\n" % (srna, start, end, strand, srna_type, feature, overlap, srnaSeq))
            else:
                srnaPredictedFile = open("%s/predicted/%s.fna" % (write_path, accession), "a")
                srnaPredictedFile.write(">%s[%s-%s,%s,%s,%s,%s]\n%s\n" % (srna, start, end, strand, srna_type, feature, overlap, srnaSeq))
                
def get_overlap_vals(subsetDat, overlaps):
    dat_len = len(subsetDat.index)
    overlapping_ids = []
    lengths = []
    start_val = 0
    end_val = 0
    for i in range(0,dat_len):
        query_val = subsetDat.iloc[i]['query_id']    
        new_start_val = min([subsetDat.iloc[i]['target_start'], subsetDat.iloc[i]['target_end']])
        new_end_val = max([subsetDat.iloc[i]['target_start'], subsetDat.iloc[i]['target_end']]) 
        if end_val > new_start_val:
            overlapping_ids.append(query_val)
            len_1 = end_val - start_val
            len_2 = new_end_val - new_start_val
            shortest_seq = min([len_1, len_2])
            overlap_start = max([start_val, new_start_val])
            overlap_end = min([end_val, new_end_val])
            overlap = (overlap_end - overlap_start)/shortest_seq
            overlaps.append(overlap)
        else:
            end_val = new_end_val
            start_val = new_start_val
            overlapping_ids = [query_val]
    return(overlaps)

def get_overlap_list(subsetDat):
    overlapping_ids = []
    overlap_list = []
    lengths = []
    start_val = 0
    end_val = 0
    shortest_seq = max(subsetDat['target_end'])
    dat_len = len(subsetDat.index)
    for i in range(0,dat_len):
        query_val = subsetDat.iloc[i]['query_id']    
        new_start_val = min([subsetDat.iloc[i]['target_start'], subsetDat.iloc[i]['target_end']])
        new_end_val = max([subsetDat.iloc[i]['target_start'], subsetDat.iloc[i]['target_end']])  
        if end_val > new_start_val:
            len_2 = new_end_val - new_start_val
            shortest_seq = min([shortest_seq, len_2])
            overlap_start = max([start_val, new_start_val])
            overlap_end = min([end_val, new_end_val])
            overlap = (overlap_end - overlap_start)/shortest_seq
            
            if overlap >= 0.5 and overlap_end - overlap_start >= 50:
                if query_val not in overlapping_ids:
                    overlapping_ids.append(query_val)
                end_val = max([end_val, new_end_val])
            else:
                overlap_list.append(overlapping_ids)
                shortest_seq = max(subsetDat['target_end'])
                end_val = new_end_val
                start_val = new_start_val
                overlapping_ids = [query_val]
        else:
            overlap_list.append(overlapping_ids)
            end_val = new_end_val
            start_val = new_start_val
            overlapping_ids = [query_val]
        if i == dat_len - 1:
            overlap_list.append(overlapping_ids)
    return(overlap_list)

def get_overlap_count(overlap_list, d):
    for l in overlap_list:
        list_len = len(l)
        if list_len == 0:
            continue
        for i in range(0,list_len - 1):
            for j in range(i+1, list_len):
                ids =[l[i], l[j]]
                ids.sort()
                current_id = "_".join(ids)
                if current_id in d:
                    d[current_id] += 1
                else:
                    d[current_id] = 1
    return(d)

def unique_set_of_overlaps(all_overlaps, ids_checked, id1, id2):
    make_new = True
    counter = 0
    if id1 in ids_checked:
        if id1 in all_overlaps:
            if id2 not in all_overlaps[id1]:
                all_overlaps[id1].append(id2)    
        else:
            counter = 0
            item_list = []
            for item in all_overlaps:
                if id1 in all_overlaps[item]:
                    item_list.append(item)
                    if id2 not in all_overlaps[item]:
                        all_overlaps[item].append(id2)
                    make_new = False
                    counter += 1
            if counter > 1:
                print(item_list[1:])
                for item in item_list[1:]:
                    for value in all_overlaps[item]:
                        if value not in all_overlaps[item_list[0]]:
                            all_overlaps[item_list[0]].append(value)
                    all_overlaps.pop(item, None)
                     
    else:
        ids_checked.append(id1)
    return(all_overlaps, ids_checked, make_new, counter)                
                
def combined_alignments(query, combined_d, ids_checked, query_matches):
    combined_ids = [query]
    ids_checked.append(query)
    max_query = query
    for i in range(0, len(query_ids)):
        ids =[query, query_ids[i]]
        ids.sort()
        current_id = "_".join(ids)
        if current_id in query_matches:
            if query_ids[i] in ids_checked:
                for key, value in combined_d.items():
                    if query_ids[i] in value:
                        max_query = key
            else:
                combined_ids.append(query_ids[i])
                ids_checked.append(query_ids[i])
                
                
    if max_query in combined_d:
        for item in combined_ids:
            if item not in combined_d[max_query]:
                combined_d[max_query].append(item)
    else:
        combined_d[max_query] = combined_ids
    return(combined_d, ids_checked)



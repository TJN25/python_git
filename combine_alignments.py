#!/Users/thomasnicholson/anaconda3/bin/python

import sys
import getopt

import pandas as pd

# from Bio import SeqIO

# import comparativesrna as srna

help_ = '''
    combine_alignments.py v 0.1 (October 2020) is a script for {}.
    Usage:
         combine_alignments.py [options] <input> <output>
    
    Options:
        -h	Display this help
        -q  Supress messages
    Input
        -i	<input> the input
        -o	<output> the output
'''


def rungetopts():
    try:
        opts, args = getopt.getopt(sys.argv[1:],
                                   "i:o:qh",
                                   ["input", "output", "quiet", "help"]
                                   )
    except getopt.GetoptError as err:
        print(err)
        print(help_)
        sys.exit(2)
    input_ = ""
    output = ""
    for o, a in opts:
        if o in ("-h", "--help"):
            print(help_)
            sys.exit()
        elif o in ("-i", "--input"):
            input_ = a
        elif o in ("-o", "--output"):
            output = a
        else:
            assert False, "unhandled option"
    if output == "":
        print("-o <output> missing. For more help_ use -h")
        sys.exit(2)
    if input_ == "":
        print("-i <input> missing. For more help_ use -h")
        sys.exit(2)
    return input_, output


def import_alignment_data(input_):
    align_dat = pd.read_csv(
        input_,
        header=None, delim_whitespace=True)

    ##split out fields that are needed for further analysis
    align_dat = align_dat.iloc[:, [1, 3]]
    align_dat.columns = ["details", "query_id"]
    align_dat[["target_contig", "coord"]] = align_dat.details.str.split("/", expand=True)
    align_dat[["target_start", "target_end"]] = align_dat.coord.str.split("-", expand=True)
    align_dat = align_dat[["query_id", "target_contig", "target_start", "target_end"]]
    align_dat["target_start"] = align_dat["target_start"].astype(str).astype(int)
    align_dat["target_end"] = align_dat["target_end"].astype(str).astype(int)
    align_dat = align_dat.sort_values(by=['target_start'])

    return align_dat


def get_overlap_vals(subset_dat, overlaps):
    dat_len = len(subset_dat.index)
    overlapping_ids = []
    start_val = 0
    end_val = 0
    for i in range(0, dat_len):
        query_val = subset_dat.iloc[i]['query_id']
        new_start_val = min([subset_dat.iloc[i]['target_start'], subset_dat.iloc[i]['target_end']])
        new_end_val = max([
            subset_dat.iloc[i]['target_start'],
            subset_dat.iloc[i]['target_end']
        ])
        if end_val > new_start_val:
            overlapping_ids.append(query_val)
            len_1 = end_val - start_val
            len_2 = new_end_val - new_start_val
            shortest_seq = min([len_1, len_2])
            overlap_start = max([start_val, new_start_val])
            overlap_end = min([end_val, new_end_val])
            overlap = (overlap_end - overlap_start) / shortest_seq
            overlaps.append(overlap)
        else:
            end_val = new_end_val
            start_val = new_start_val
            overlapping_ids = [query_val]
    return overlaps


def get_overlap_list(subset_dat_):
    overlapping_ids = []
    overlap_list = []
    start_val = 0
    end_val = 0
    shortest_seq = max(subset_dat_['target_end'])
    dat_len = len(subset_dat_.index)
    for i in range(0, dat_len):
        query_val = subset_dat_.iloc[i]['query_id']
        new_start_val = min([subset_dat_.iloc[i]['target_start'], subset_dat_.iloc[i]['target_end']])
        new_end_val = max([subset_dat_.iloc[i]['target_start'], subset_dat_.iloc[i]['target_end']])
        if end_val > new_start_val:
            len_2 = new_end_val - new_start_val
            shortest_seq = min([shortest_seq, len_2])
            overlap_start = max([start_val, new_start_val])
            overlap_end = min([end_val, new_end_val])
            overlap = (overlap_end - overlap_start) / shortest_seq

            if overlap >= 0.5 and overlap_end - overlap_start >= 50:
                overlapping_ids.append(query_val)
                end_val = max([end_val, new_end_val])
            else:
                overlap_list.append(overlapping_ids)
                shortest_seq = max(subset_dat_['target_end'])
                end_val = new_end_val
                start_val = new_start_val
                overlapping_ids = [query_val]
        else:
            overlap_list.append(overlapping_ids)
            end_val = new_end_val
            start_val = new_start_val
            overlapping_ids = [query_val]
    return overlap_list


def get_overlap_count(overlap_list, d):
    for item in overlap_list:
        list_len = len(item)
        if list_len == 0:
            continue
        for i in range(0,list_len - 1):
            for j in range(i+1, list_len):
                ids =[item[i], item[j]]
                ids.sort()
                current_id = "_".join(ids)
                if current_id in d:
                    d[current_id] += 1
                else:
                    d[current_id] = 1
    return(d)


def combined_alignments(query, combined_d, ids_checked, query_matches):
    combined_ids = [query]
    ids_checked.append(query)
    max_query = query
    for i in range(0, len(query_ids)):
        ids = [query, query_ids[i]]
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
    return (combined_d, ids_checked)


def main():
    infile, outfile = rungetopts()

    ## Import data and get list of queries and targets
    align_dat = import_alignment_data(input_=infile)
    target_contigs = align_dat.target_contig.unique()
    query_ids = align_dat.query_id.unique()

    combined_ids_count = {}

    for contig in target_contigs:
        print(contig)
        subset_dat = align_dat.loc[align_dat['target_contig'] == contig]
        overlap_list = get_overlap_list(subset_dat_=subset_dat)
        combined_ids_count = get_overlap_count(overlap_list=overlap_list, d=combined_ids_count)
    print(combined_ids_count)

    target_counts = {}
    for query in query_ids:
        print(query)
        query_dat = align_dat.loc[align_dat['query_id'] == query]
        query_len = len(query_dat.index)
        target_counts[query] = query_len

    query_matches = {}
    overlap_percentages = []
    for i in range(0, len(query_ids) - 1):
        for j in range(i + 1, len(query_ids)):
            ids = [query_ids[i], query_ids[j]]
            ids.sort()
            current_id = "_".join(ids)
            count = min([target_counts[query_ids[i]], target_counts[query_ids[j]]])
            if current_id in combined_ids_count:
                if combined_ids_count[current_id] / count < 1:
                    overlap_percentages.append(combined_ids_count[current_id] / count)
                else:
                    overlap_percentages.append(1)
                if combined_ids_count[current_id] / count > 0.35:
                    query_matches[current_id] = combined_ids_count[current_id] / count

    combined_d = {}
    ids_checked = []
    for query in query_ids:
        if query not in ids_checked:
            combined_d, ids_checked = combined_alignments(query=query, combined_d=combined_d, ids_checked=ids_checked,
                                                          query_matches=query_matches)

    print(combined_d)

if __name__ == "__main__":
    main()

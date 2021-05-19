import os, time, datetime
import re
import json
import requests
import argparse
import logging
logging.basicConfig(level=logging.INFO)

def jprint(obj):
    # create a formatted string of the Python JSON object
    text = json.dumps(obj, sort_keys=True, indent=4)
    print(text)
    open(file="test_out.json", mode="w").write(text)

def get_latest_pango_lineages_list():
    response = requests.get("https://raw.githubusercontent.com/cov-lineages/pango-designation/master/lineage_notes.txt")
    pango_lineages_list=[]
    for row in response.iter_lines():
        lineage = row.decode("utf8").strip().split("\t")[0]
        if '*' not in lineage:
            pango_lineages_list.append(lineage)
    return pango_lineages_list[1:] #the first element is title 'Lineage' so we return from index 1 and not 0

def get_mutation_info_from_api():
    pango_lineages_list = get_latest_pango_lineages_list()
    dict_of_mutation_hits_by_lineage = {}

    n_pango_lineages = len(pango_lineages_list)
    for i, pango_lineage in enumerate(pango_lineages_list):
        logging.info("Caching mutation data for {} ({}/{})".format(pango_lineage, i + 1, n_pango_lineages))
        response = requests.get(
            "https://api.outbreak.info/genomics/lineage-mutations?pangolin_lineage=" + pango_lineage + "&frequency=0")
        dict_of_mutation_hits_by_lineage[pango_lineage] = response.json()["results"]

    with open('lineage_muation_cache.txt', 'w') as file:
        file.write(json.dumps(
            dict_of_mutation_hits_by_lineage))  # use to load cache json.load(open(file="file.txt", mode="r"))

def tsv_results_output_render(results_single_dict, results_multi_dict):

    day_stamp = datetime.date.today()
    out_tsv_filename="result_mutations_in_pangolineages_{}-{}-{}.tsv".format(
                                                       day_stamp.year,
                                                       day_stamp.month,
                                                       day_stamp.day)

    #Header
    with open(file=out_tsv_filename, mode="w") as fp_out:
        fp_out.write("Mutation\t#PangoLineages\tPangoNames(prevalence in lineage)\n")


    for query_mutation in results_single_dict.keys():
        with open(file=out_tsv_filename, mode="a") as fp_out:
            lineages_prevalence_sorted_list=sorted([(key, results_single_dict[query_mutation][key]['prevalence'])
                                                    for key in results_single_dict[query_mutation].keys()],
                                                   key=lambda x: x[1],
                                                   reverse=True)

            if not lineages_prevalence_sorted_list:
                lineages_prevalence_sorted_list=[("NA","NA")]


            fp_out.write("{}\t{}\t{}\n".format(query_mutation,
                             len(results_single_dict[query_mutation].keys()),
                             ",".join(["{}({})".format(item[0], item[1]) for item in lineages_prevalence_sorted_list])
                                             )
                     )

    for query_multi_mutation in results_multi_dict.keys():
        #print(results_multi_dict[query_multi_mutation])
        with open(file=out_tsv_filename, mode="a") as fp_out:
            fp_out.write("{}\t{}\t{}\n".format(query_multi_mutation,
                                               len(results_multi_dict[query_multi_mutation].keys()),
                                               ",".join(results_multi_dict[query_multi_mutation].keys()))

                         )
    logging.info("Results written to file: {}".format(out_tsv_filename))

def csv_list(string):
   return string.split(',')

def cli_args():
    parser = argparse.ArgumentParser(
        "{} uses https://api.outbreak.info/ to get information on query mutations supported by GISAID database data\n".format(os.path.basename(__file__)))

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-m', '--mutation_names', type = csv_list,
                        required=False, help="Mutation name to query (e.g. s:del69/70, s:d614g)"
                        )
    group.add_argument('-f', '--file_mutation_list', required=False,
                       help="Single column input file with mutation lists")

    return parser.parse_args()

def multi_mutation_query(query_multi_mutations,pango_lineages_list, mutations_lineages_dict):
    results_dict={}
    for query_multi_mutation in query_multi_mutations:

        results_dict[query_multi_mutation] = {}
        for pango_lineage in pango_lineages_list:
            results_dict[query_multi_mutation][pango_lineage] = {}
            query_multi_mutation_presence = {i: False for i in query_multi_mutation.split("+")}
            mutations_dict = mutations_lineages_dict[pango_lineage]

            for mutation_meta_dict in mutations_dict:
                if mutation_meta_dict["mutation"] in query_multi_mutation_presence.keys():
                    #print("{}:{}".format(pango_lineage,mutation_meta_dict["mutation"]))
                    query_multi_mutation_presence[mutation_meta_dict["mutation"]]=True
                    results_dict[query_multi_mutation][pango_lineage][mutation_meta_dict['mutation']]={}
                    results_dict[query_multi_mutation][pango_lineage][mutation_meta_dict['mutation']]['prevalence']=mutation_meta_dict['prevalence']
                    results_dict[query_multi_mutation][pango_lineage][mutation_meta_dict['mutation']]['mutation_count'] = mutation_meta_dict['mutation_count']
                    results_dict[query_multi_mutation][pango_lineage][mutation_meta_dict['mutation']]['lineage_count'] = mutation_meta_dict['lineage_count']
            if not all(query_multi_mutation_presence.values()):
                #print("Both mutations are NOT present in {}. Removing from results".format(pango_lineage))
                del results_dict[query_multi_mutation][pango_lineage]

    return results_dict

def parse_input_text_file(filepath):
    with open(filepath, mode="r") as fp_in:
        return fp_in.readlines()

def single_mutation_query(query_single_mutations, pango_lineages_list, mutations_lineages_dict):
    results_dict = {}
    for query_mutation in query_single_mutations:
        results_dict[query_mutation] = {}
        for pango_lineage in pango_lineages_list:
            mutations_dict = mutations_lineages_dict[pango_lineage]
            results_dict[query_mutation][pango_lineage] = {}
            for mutation_meta_dict in mutations_dict:
                if mutation_meta_dict["mutation"] == query_mutation:
                    results_dict[query_mutation][pango_lineage]['prevalence'] = mutation_meta_dict['prevalence']
                    results_dict[query_mutation][pango_lineage]['mutation_count'] = mutation_meta_dict['mutation_count']
                    results_dict[query_mutation][pango_lineage]['lineage_count'] = mutation_meta_dict['lineage_count']

            if not results_dict[query_mutation][pango_lineage].keys():
                del results_dict[query_mutation][pango_lineage]
    return results_dict


def main():
    args=cli_args()

    if not os.path.exists('lineage_muation_cache.txt') or (time.time()-os.stat('lineage_muation_cache.txt').st_mtime) > 3600*24*5:
        get_mutation_info_from_api()

    if args.mutation_names:
        query_mutations=[mutation.lower() for mutation in args.mutation_names]
    else:
        query_mutations=[i.rstrip() for i in parse_input_text_file(args.file_mutation_list)]



    bool_multi_mutations = ["+" in i for i in query_mutations]
    query_single_mutations = [i[0] for i in zip(query_mutations,bool_multi_mutations) if i[1] == False]
    query_multi_mutations = [i[0] for i in zip(query_mutations, bool_multi_mutations) if i[1] == True]

    mutations_lineages_dict = json.load(open(file="lineage_muation_cache.txt", mode="r"))
    pango_lineages_list=mutations_lineages_dict.keys()

    results_single_dict = single_mutation_query(query_single_mutations,
                                        pango_lineages_list,
                                         mutations_lineages_dict)
    results_multi_dict = multi_mutation_query(query_multi_mutations,pango_lineages_list, mutations_lineages_dict)


    tsv_results_output_render(results_single_dict,results_multi_dict)

if __name__ == '__main__':
    main()
    print("Done")



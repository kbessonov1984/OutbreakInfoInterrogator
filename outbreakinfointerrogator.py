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

def tsv_results_output_render(results_dict, mutation_name):
    mutation_name_norm = re.sub('/', '-', mutation_name)
    day_stamp = datetime.date.today()
    out_tsv_filename="result_mutation_{}_{}-{}-{}.tsv".format(mutation_name_norm,
                                                       day_stamp.year,
                                                       day_stamp.month,
                                                       day_stamp.day)

    lineages_prevalence_sorted_list=sorted([(key, results_dict[key]['prevalence']) for key in results_dict.keys()],
                                           key=lambda x: x[1],
                                           reverse=True)

    with open(file=out_tsv_filename, mode="w") as fp_out:
        fp_out.write("Mutation\t#PangoLineages\tPangoNames(prevalence in lineage)\n")
        fp_out.write("{}\t{}\t{}".format(mutation_name,
                         len(results_dict.keys()),
                         ",".join(["{}({})".format(item[0], item[1]) for item in lineages_prevalence_sorted_list])
                                         )
                     )
    logging.info("Results written to file: {}".format(out_tsv_filename))

def cli_args():
    parser = argparse.ArgumentParser(
        "{} uses https://api.outbreak.info/ to get information on query mutations supported by GISAID database data\n".format(os.path.basename(__file__)))

    parser.add_argument('-m', '--mutation_name',
                        required=True, help="Mutation name to query (e.g. s:del69/70, s:d614g, orf1a:del17)"
                        )
    return parser.parse_args()

def main():
    args=cli_args()

    if not os.path.exists('lineage_muation_cache.txt') or (time.time()-os.stat('lineage_muation_cache.txt').st_mtime) > 86400:
        get_mutation_info_from_api()

    query_mutation=args.mutation_name
    #query_mutation="s:d614g"
    #query_mutation="orf1a:del17"
    #query_mutation="s:del69/70"
    mutations_lineages_dict = json.load(open(file="lineage_muation_cache.txt", mode="r"))
    pango_lineages_list=mutations_lineages_dict.keys()

    results_dict={}
    for pango_lineage in pango_lineages_list:
        mutations_dict = mutations_lineages_dict[pango_lineage]
        results_dict[pango_lineage]={}
        for muation_meta in mutations_dict:
            if muation_meta["mutation"] == query_mutation:
                results_dict[pango_lineage]['prevalence'] = muation_meta['prevalence']
                results_dict[pango_lineage]['mutation_count'] = muation_meta['mutation_count']
                results_dict[pango_lineage]['lineage_count'] = muation_meta['lineage_count']

        if not results_dict[pango_lineage].keys():
            del results_dict[pango_lineage]




    tsv_results_output_render(results_dict, query_mutation)

if __name__ == '__main__':
    main()
    print("Done")

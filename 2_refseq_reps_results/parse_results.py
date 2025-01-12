## Given a directory of antismash results, extract NRPS and NRP-metallophore clusters

import os
from sys import argv
import re
import itertools
import csv

from Bio import SeqIO
import pandas as pd

from tqdm import tqdm

TARGET_HMMS = ["VibH_like", "Cy_tandem",
               "IBH_Asp", "IBH_His", "TBH_Asp", "SBH_Asp",
               "CyanoBH_Asp1", "CyanoBH_Asp2",
               "IPL", "SalSyn", "EntA", "EntC",
               "GrbD", "GrbE",
               "FbnL", "FbnM",
               "PvdO", "PvdP",
               "Orn_monoox", "Lys_monoox", "VbsL",
               "KtzT", "MetRS-like",
               "FecCD", "TonB_dep_Rec", "Peripla_BP_2"
               ]


def parse_domains(gbpath):
    try:
        og_region = re.findall("region\d\d\d", gbpath)[0]
        assembly = re.findall("GCF_[0-9.]+", gbpath)[0]
    except:
        og_region = ""
        assembly = gbpath

    record = SeqIO.read(gbpath, "genbank")
    hits = set()
    hitscores = []
    original_id = record.annotations["structured_comment"]["antiSMASH-Data"].get("Original ID")
    if not original_id:
        record_id = record.id
    else:
        record_id = original_id.split(" ")[0]

    start = record.annotations["structured_comment"]["antiSMASH-Data"]['Orig. start']
    end = record.annotations["structured_comment"]["antiSMASH-Data"]['Orig. end']
    taxonomy = record.annotations["taxonomy"]
    organism = record.annotations["organism"]
    contig_edge = "MISSING"
    for feat in record.features:
        if feat.type == "region":
            contig_edge = feat.qualifiers.get("contig_edge")[0]
            products = feat.qualifiers.get("product")
            nrps = str("NRPS" in products)
            nrpslike = str("NRPS-like" in products)
            nrpsid = str("NRP-siderophore" in products)
            # region = feat.qualifiers.get("region_number")
        accession = feat.qualifiers.get("protein_id")
        if not accession:
            accession = feat.qualifiers.get("locus_tag")
        if not accession:
            accession = feat.qualifiers.get("ID")
        domains = feat.qualifiers.get("sec_met_domain")
        if domains:
            hmms = [dom.split(" ")[0] for dom in domains]
            bitscores = [re.findall(r'(?<=bitscore: )[0-9.]+', dom)[0] for dom in domains]
            hitscores.extend(zip(hmms, bitscores, itertools.cycle(accession)))
            for h in hmms:
                if h in TARGET_HMMS:
                    hits.add(h)
    # Assembly, Accession, Region, NRPS, NRPS-like, NRP-siderophore, Span, ContigEdge, Taxonomy, Species, HMMs, HMM scores
    return [assembly, record_id, og_region, nrps, nrpslike, nrpsid, "..".join((start, end)), contig_edge, "_".join(taxonomy), organism,
            ",".join(hits), hitscores]
    # return {"Accession": record.id, "Taxonomy": record.annotations.get("taxonomy"), "Species": record.annotations.get("organism"), "HMMs": " ".join(hits)}


# Given a list of lists of tuples, clean it up into a list of dictionaries
def clean_scores(scores_list):
    top_hits = []
    for scores in scores_list:
        hit_dict = {}
        for hmm, score, acc in scores:  # accession isn't used
            if hmm not in TARGET_HMMS:
                continue
            # New HMM
            if hmm not in hit_dict:
                hit_dict[hmm] = score
            # Existing HMM goes to winner
            else:
                if float(hit_dict[hmm]) < float(score):
                    hit_dict[hmm] = score
        top_hits.append(hit_dict)
    return top_hits


def gene_tab(record_infos, outpath):
    all_hmms = set()
    for record in record_infos:
        all_hmms.update(record[14].split(","))
    all_hmms = list((h for h in all_hmms if h != ""))

    row_list = []
    for record in record_infos:
        r_ass = record[0]
        r_acc = record[1]
        r_region = record[2]
        r_taxon = record[5]
        r_species = record[6]
        scores = record[-1]

        cluster_domains = {}
        for hmm, score, acc in scores:
            if hmm not in TARGET_HMMS:
                continue
            if acc not in cluster_domains:
                cluster_domains[acc] = {}
            cluster_domains[acc][hmm] = score

        for key, val in cluster_domains.items():
            core_info = [r_ass, r_acc, r_region, r_taxon, r_species, key]
            domain_info = [val.get(dom, "0") for dom in all_hmms]
            row_list.append(core_info + domain_info)

    with open(outpath, "w") as outf:
        handle = csv.writer(outf, delimiter="\t")
        handle.writerow(["Assembly", "Accession", "Region", "Taxonomy", "Species", "Protein"] + all_hmms)
        for row in row_list:
            handle.writerow(row)


def parse_results(asdir):
    record_infos = []
    for genome in tqdm(os.listdir(asdir)):
        if genome.endswith(".tsv"):
            continue

        genome_dir = os.path.join(asdir, genome)
        regions = (r for r in os.listdir(genome_dir) if "region" in r)

        for region_path in regions:
            # Do basic parsing to find if it's relevant
            relevant = False
            region = False
            with open(os.path.join(genome_dir, region_path), 'r') as record:
                try:
                    for line in record:
                        if region == True:
                            if line.strip().endswith(
                                    ('/product="NRP-siderophore"', '/product="NRPS"', '/product="NRPS-like"')):
                                relevant = True
                                break
                            if ".." in line:  # Next feature
                                break
                        if line.startswith("     region"):
                            region = True
                except UnicodeDecodeError:
                    print(region_path)
                    raise
            if relevant:
                record_info = parse_domains(os.path.join(genome_dir, region_path))
                record_infos.append(record_info)


    # Convert to dataframe (except last column)
    df = pd.DataFrame((r[:-1] for r in record_infos),
                      columns = ("Assembly", "Accesion", "Region", "NRPS", "NRPS-like", "NRP-siderophore", "Span", "ContigEdge", "Taxonomy", "Species", "HMMs"))


    print("Preparing score table")
    scores = [info[-1] for info in record_infos]
    cleaned = clean_scores(scores)
    score_df = pd.DataFrame.from_dict(cleaned)
    df_combo2 = pd.concat([df.iloc[:,0:-1], score_df], axis = 1)

    print("Writing score table")
    with open(os.path.join(asdir, "scores.tsv"), 'w') as outf:
        df_combo2.to_csv(outf, sep = "\t")

    print("Writing gene table")
    gene_tab(record_infos, os.path.join(asdir, "genes.tsv"))

if __name__ == "__main__":
    asdir = argv[1]
    parse_results(asdir)

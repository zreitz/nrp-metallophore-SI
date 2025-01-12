# Find the intersection between complete RefSeq representative genomes and the GTDB tree
from ete3 import Tree
import csv

# RefSeq representative genomes
with open("../2_refseq_reps_results/assembly_summary.txt", 'r') as inf:
    summary = list(csv.reader(inf, delimiter = "\t"))

complete = set((l[0][4:13] for l in summary if l[11] == "Complete Genome"))

# This tree came from https://data.gtdb.ecogenomic.org/releases/latest/
t = Tree('bac120.tree', quoted_node_names=True, format=1)

# Rename leaves to just the 9 digits
leaves = t.get_leaves()
for i in range(len(leaves)):
    leaves[i].name = leaves[i].name[7:16]
in_tree = set((l.name for l in leaves))

# Keep only leaves that were scanned (and in the tree)
t.prune(complete.intersection(in_tree))
print(len(t.get_leaves()))

t.write(format=0, outfile = "pruned.tree")

taxa_strs = [leaf.name for leaf in t.get_leaves()]

with open("tree_taxa.txt", 'w') as outf:
    for leaf in taxa_strs:
        outf.write(leaf + "\n")

with open("assembly_summary.txt", 'r') as inf, open("tree_assembly_summary.txt", 'w') as outf:
    for line in inf.readlines():
        if line[4:13] in taxa_strs:
            outf.write(f"{line[4:13]}\t{line}")
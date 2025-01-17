
#!/usr/bin/python3.8
# Created by Bita

##reads a tree in newick format from a text file, loads it as ete3 tree
##reads a list of nodes from files, and prunes the tree to only keep the leaves from the list

import re
import pandas as pd
import os
import ete3
import argparse


# GTDB_tree_string = open('bac120.tree', 'r').read()
# print(GTDB_tree)


def main(input_tree, labels_file, output):


    ete_tree = ete3.PhyloTree(input_tree, format=1, quoted_node_names=True)

#    ete_tree.write(format=1, outfile="GTDB_tree_ete3.nwk")
    # print(ete_tree)

#    collapse_tree = ete_tree.collapse_lineage_specific_expansions()
#    print(collapse_tree)
#    breakpoint()
##traversing through all nodes (both leaves and internal nodes)
# for node in ete_tree.traverse("postorder"):
# Do some analysis on node
#  print(node.name, node.is_leaf())


#    user_prune_file = input("Enter the file name with the prune list: ")
    prune_file = open(labels_file, 'r')
    prune_list = prune_file.read().splitlines()
    prune_file.close()



    # print(type(prune_list))
    print("pruned tree is: \n")
    ete_tree.prune(prune_list, preserve_branch_length=False)
    print(ete_tree)


    # Write the pruned tree to the output file
    ete_tree.write(format=1, outfile=output)
    # Write the pruned tree to the output file
# Parsing arguments & Execution #

if __name__ == "__main__":
    # Parsing arguments
    parser = argparse.ArgumentParser(description='Takes reference tree and list of nodes and Prunes the tree')
    parser.add_argument('-t', '--ref_tree', help='Reference Tree', required=True)
    parser.add_argument('-l', '--labels', help='File with labels', required=True)
    parser.add_argument('-o', '--output', help='Output tree newick file', required=True)
    args = parser.parse_args()

    # Execution
    main(args.ref_tree, args.labels, args.output)

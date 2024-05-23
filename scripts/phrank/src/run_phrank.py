import sys
import os
from collections import defaultdict
from phrank import *

DAGFILE = f"/run/data_dependencies/phrank/{sys.argv[3]}/child_to_parent.txt"
DISEASEANNOTATIONS = f"/run/data_dependencies/phrank/{sys.argv[3]}/disease_to_pheno.txt"
GENEANNOTATIONS = f"/run/data_dependencies/phrank/{sys.argv[3]}/gene_to_phenotype.txt"
DISEASEGENE = f"/run/data_dependencies/phrank/{sys.argv[3]}/disease_to_gene.txt"

#/houston_10t/cole/MARRVEL/udn/phrank/data/knowledgebase/gene_to_symbol.txt
#/houston_10t/cole/MARRVEL/udn/phrank/data/knowledgebase/hpo_to_name.txt

def load_list(listfile):
  returnList = []
  for line in open(listfile):
    returnList.append(line.strip())
  return returnList

def load_set(listfile):
  returnSet = set()
  for line in open(listfile):
    returnSet.add(line.strip().split("\t")[0])
  return returnSet

def main(genelist, phenotypelist):
  p = Phrank(DAGFILE, DISEASEANNOTATIONS, DISEASEGENE, GENEANNOTATIONS)
  genes = load_set(genelist)
  phenos = load_set(phenotypelist)
  ranking = p.rank_genes(genes, phenos)
  for item in ranking: print(str(item[1]) + "\t" + str(item[0]))

if __name__ == "__main__": main(sys.argv[1], sys.argv[2])


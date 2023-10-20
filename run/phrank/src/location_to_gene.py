import sys
import os
from collections import defaultdict


if sys.argv[2] == 'hg19':
    GENE_LOCS = "/run/data_dependencies/phrank/hg19/grch37_symbol_to_location.txt"
else:
    GENE_LOCS = "/run/data_dependencies/phrank/hg38/grch38_symbol_to_location.txt"

def load_loc_maps(locmap):
  by_start = defaultdict(list)
  by_end = defaultdict(list)
  by_symbol = defaultdict(lambda: defaultdict(lambda: None))
  for line in open(locmap):
    lineData = line.strip().split("\t")
    symbol = lineData[0]
    if "." in symbol: continue
    chrom = lineData[1]
    start = int(lineData[2])
    end = int(lineData[3])
    if not by_symbol[chrom][symbol]: by_symbol[chrom][symbol] = [start, end]
    else:
        if by_symbol[chrom][symbol][0] > start: by_symbol[chrom][symbol][0] = start
        if by_symbol[chrom][symbol][1] < end: by_symbol[chrom][symbol][1] = end
    by_start[chrom].append((start, end, symbol))
    by_end[chrom].append((end, start, symbol))
  for chrom in by_start.keys():
    by_start[chrom].sort()
    by_end[chrom].sort()
  return by_start, by_end, by_symbol

def binary_search(loc_to_symbol, start, end):
  minimum = 0 #Inclusive
  maximum = len(loc_to_symbol) #Exclusive
  returnSet = set()
  index = None
  while minimum < maximum:
    index = (minimum + maximum) // 2
    location = loc_to_symbol[index]
    if location[0] < start:
      minimum = index + 1
      continue
    if location[0] > end:
      maximum = index
      continue
    break
  if index is None: return returnSet
  returnSet.add(loc_to_symbol[index][2])
  i = index - 1
  while i >= 0:
    location = loc_to_symbol[i]
    i -= 1
    if location[0] < start: break
    returnSet.add(location[2])
  i = index + 1
  while i < len(loc_to_symbol):
    location = loc_to_symbol[i]
    i += 1
    if location[0] > end: break
    returnSet.add(location[2])
  return returnSet

def find_overlapping_genes(by_symbol, by_start, by_end):
  gene_overlaps = defaultdict(set)
  for chrom in by_symbol.keys():
    for symbol in by_symbol[chrom].keys():
      start, end = by_symbol[chrom][symbol]
      start_overlaps = binary_search(by_start[chrom], start, end)
      end_overlaps = binary_search(by_end[chrom], start, end)
      gene_overlaps[symbol] = (start_overlaps | end_overlaps) - set([symbol])
  return gene_overlaps

def main(variants):
  by_start, by_end, by_symbol = load_loc_maps(GENE_LOCS)
  #atxn1start = by_symbol["6"]["ENSG00000124788"][0]
  #atxn1end = by_symbol["6"]["ENSG00000124788"][1]
  for line in open(variants):
    lineData = line.strip().split(":")
    chrom = lineData[0]
    pos = int(lineData[1])
    start_overlaps = binary_search(by_start[chrom], pos, pos)
    end_overlaps = binary_search(by_end[chrom], pos, pos)
    genes = start_overlaps | end_overlaps
    for gene in genes: print(line.strip() + "\t" + gene)

if __name__ == "__main__":
  main(sys.argv[1])


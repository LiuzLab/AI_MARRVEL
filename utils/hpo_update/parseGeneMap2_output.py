#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# This script parses the genemap2.txt file from OMIM and outputs
# all fields along with processed phenotype information.
# Usage: cat genemap2.txt | ./parseGeneMap2.py > parsed.txt
#

#originially from: https://github.com/OMIM-org/genemap2-parser/tree/master

# Imports
import sys
import re

# Define the expected number of fields in genemap2.txt
EXPECTED_NUM_FIELDS = 14

# Print header for the output file (TSV format)
header = [
    "Chromosome",
    "Genomic_Position_Start",
    "Genomic_Position_End",
    "Cyto_Location",
    "Computed_Cyto_Location",
    "MIM_Number",
    "Gene_Symbols",
    "Gene_Name",
    "Approved_Gene_Symbol",
    "Entrez_Gene_ID",
    "Ensembl_Gene_ID",
    "Comments",
    "Mouse",
    "Phenotype_Name",
    "Phenotype_MIM_Number",
    "Mapping_Key",
    "Inheritance"
]
print("\t".join(header))

# Read from stdin
for line_number, line in enumerate(sys.stdin, 1):
    # Skip comments
    if line.startswith('#'):
        continue

    # Strip trailing new line and carriage return
    line = line.strip('\n').strip('\r')

    # Split the line into fields based on tab
    valueList = line.split('\t')

    # Validate the number of fields
    if len(valueList) < EXPECTED_NUM_FIELDS:
        sys.stderr.write(f"Warning: Line {line_number} has fewer fields ({len(valueList)}) than expected ({EXPECTED_NUM_FIELDS}). Skipping.\n")
        continue  # Skip malformed lines

    # Assign fields to variables
    chromosome = valueList[0]
    genomicPositionStart = valueList[1]
    genomicPositionEnd = valueList[2]
    cytoLocation = valueList[3]
    computedCytoLocation = valueList[4]
    mimNumber = valueList[5]
    geneSymbols = valueList[6]
    geneName = valueList[7]
    approvedGeneSymbol = valueList[8]
    entrezGeneID = valueList[9]
    ensemblGeneID = valueList[10]
    comments = valueList[11]
    phenotypeString = valueList[12]
    mouse = valueList[13]

    # Skip empty phenotypes
    if not phenotypeString:
        continue

    # Parse the phenotypes
    phenotypes = phenotypeString.split(';')
    for phenotype in phenotypes:
        # Clean the phenotype string
        phenotype = phenotype.strip()

        # Initialize phenotype-related variables
        phenotype_name = ""
        phenotype_mim_number = ""
        phenotype_mapping_key = ""
        inheritances = ""

        # Attempt to match the long phenotype pattern
        matcher = re.match(r'^(.*),\s(\d{6})\s\((\d)\)(?:,\s*(.*))?$', phenotype)
        if matcher:
            phenotype_name = matcher.group(1).strip()
            phenotype_mim_number = matcher.group(2).strip()
            phenotype_mapping_key = matcher.group(3).strip()
            inheritances_raw = matcher.group(4)

            if inheritances_raw:
                inheritances = ','.join([inh.strip() for inh in inheritances_raw.split(',')])

        else:
            # Attempt to match the short phenotype pattern
            matcher = re.match(r'^(.*)\((\d)\)(?:,\s*(.*))?$', phenotype)
            if matcher:
                phenotype_name = matcher.group(1).strip()
                # phenotype_mim_number remains empty for short phenotypes
                phenotype_mapping_key = matcher.group(2).strip()
                inheritances_raw = matcher.group(3)

                if inheritances_raw:
                    inheritances = ','.join([inh.strip() for inh in inheritances_raw.split(',')])

        # Prepare the output row
        output_fields = [
            chromosome,
            genomicPositionStart,
            genomicPositionEnd,
            cytoLocation,
            computedCytoLocation,
            mimNumber,
            geneSymbols,
            geneName,
            approvedGeneSymbol,
            entrezGeneID,
            ensemblGeneID,
            comments,
            mouse,
            phenotype_name,
            phenotype_mim_number,
            phenotype_mapping_key,
            inheritances
        ]

        # Replace any None or missing inheritances with empty string
        output_fields = [field if field else "" for field in output_fields]

        # Join fields with tabs and print
        print("\t".join(output_fields))

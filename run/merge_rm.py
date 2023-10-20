import pandas as pd
import sys

###merge score files with final rank matrix 

#read rank matrix and score file
rm=pd.read_csv("/out/final_matrix/"+sys.argv[1]+".csv")
annot=pd.read_csv("/out/scores/"+sys.argv[1]+".txt.gz", sep="\t", compression="gzip")

#annot=pd.read_csv("/local/LucianLi/rami-test/r1_"+sys.argv[1]+"_scores.txt", sep="\t", compression="gzip")

#get columns from score file
annot=annot[["varId","varId_dash", "geneSymbol", "geneEnsId", "rsId", "HGVSc", "HGVSp", "IMPACT", "Consequence", "phenoList", "phenoInhList", "clin_code", "clinvarCondition"]]

#get original coordinates
annot["varId"]=annot["varId"].apply(lambda x: x.split("_E")[0])

#rename for final output
annot=annot.rename(columns={"varId":"origId", "IMPACT": "IMPACT_text", "clin_code": "clinvarSignDesc"})

#merge
test=rm.merge(annot, right_on="varId_dash", left_on="Unnamed: 0", how="left")

test.to_csv("/out/final_matrix_expanded/"+sys.argv[1]+".expanded.csv.gz", index=False, compression="gzip")


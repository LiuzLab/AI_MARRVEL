.. _output:

*********************
AIM Output Details
*********************


Output Files
========================
Following the execution of AIM, you will receive five output files, which are saved to the designated output directory.

.. list-table:: AIM output files
   :header-rows: 1
   :widths: 40 60
   
   *  -  File Name
      -  Explanation
   *  -  \*_default_predictions.csv
      -  Prediction results from AIM Default mode.
   *  -  \*_recessive_predictions.csv
      -  Prediction results from AIM Recessive mode. Each row is a pair of variants.
   *  -  \*_nd_predictions.csv
      -  Prediction results from AIM Novel Disease Gene mode.
   *  -  \*_nd_recessive_predictions.csv
      -  Prediction results from AIM Novel Disease Gene + Recessive mode.
   *  -  \*_integrated.csv
      -  Combining all 4 prediction results.

Output Explanation
======================
Each of the output files mentioned above, contains:

- Variant ID
- Annotation information (features)
- Prediction results: score, rank, and confidence
- Extra information about the variant (only in \*_integrated.csv)




Prediction Results
-----------------------

``predict``
   Prediction score from AIM

``ranking``
   Maximum ranking of all variants by prediction scores

``confidence``
   Confidence score. It's the z-score of AIM score in our inner cohort.

``confidence level``
   High / Medium / Low. We set confidence score cutoff based on precision and recall.




Extra information
-----------------------
``geneSymbol``, ``geneEnsId``, ``rsId``
   Variant and Gene ID

``HGVSc``
   Human Genome Variation Society's nomenclature for describing variants at the cDNA level.

``HGVSp``
   Human Genome Variation Society's nomenclature for describing variants at protein level.




Variant Annotation / AIM Features
-------------------------------------

AIM takes information from various database and perform feature engineering, to make the final prediction.
Here, we provide both raw and engineered features as part of the output.

In making the decision, AIM takes various types of information into consideration.

**Conservation**

``GERPpp_RS``
   GERP++ RS score, the larger the score, the more conserved the site. Scores range from -12.3 to 6.17

``LRT_Omega``
   Estimated nonsynonymous-to-synonymous-rate ratio (Omega, reported by LRT)

``LRT_score``
   The original LRT two-sided p-value (LRTori), ranges from 0 to 1.

``phyloP100way_vertebrate``
   phyloP (phylogenetic p-values) conservation score based on the multiple alignments of 100 vertebrate genomes (including human). The larger the score, the more conserved the site. Scores range from -20.0 to 10.003 in dbNSFP.

**Constrain**

``hom``
   Number of homozygotes variant in gnomAD 

``decipherVarFound``
   0/1. Whether the variant is found in the a deletion of the DECIPHER control database.

``dgvVarFound``
   0/1. Whether the variant is found in a deletion of the DGV database.

``conservationScoreDGV``
   1/3 (Low/High). If DGV subtype is Loss or Deletion, score will be 1. Otherwise 3.

``gnomadGeneOELof``
   observed/expected ratio of loss-of-function variants in gnomAD database.

``gnomadGeneOELofUpper``
   The upper bound of the confidence interval for OE LoF.

``conservationScoreOELof``
   1/2 (Low/High). If ``gnomadGeneOELofUpper < 0.35``, score is 1; otherwise 2.

``gnomadGenePLI``
   pLI score stands for the "probability of being loss-of-function intolerant" in gnomAD.

``gnomadGeneZscore``
   The gene z-score in gnomAD is related to missense variants and reflects how many standard deviations the observed count of missense variants is from the expected count.
   This metric can help identify genes that are under selective pressure and may be related to diseases.

``conservationScoreGnomad``
   1/2 (Low/High). If both gnomadAF and gnomadAFg are less than 0.01, score is high; otherwise low.



**Disease Database**

``CLASS``
   CLASS from HGMD
      DM: disease-causing mutation;
      DM? Likely disease-causing, but with questionable pathogenicity

``clinVarGeneFound``
   0/1, whether or not variant gene is found in ClinVar.

``clinVarVarFound``
   0/1, whether or not variant itself is found in ClinVar.

``curationScoreClinVar``
   1/2/3 (Low/Medium/High), curated using ClinVar significance description.

``isB/LB``
   0/1. Whether ClinVar significance description contains benign and no conflicting interpretation.

``isP/LP``
   Float ranging 0-1. Among all descriptions in ClinVar about this variant, proportion of pathogenic ones.

``clinvarNumB``
   Proportion of benign variants in the variant gene.

``clinvarNumLB``
   Proportion of likely benign + benign variants in the variant gene.

``clinvarNumLP``
   Proportion of likely pathogenic + pathogenic variants in the variant gene.

``clinvarNumP``
   Proportion of pathogenic variants in the variant gene.

``hgmdGeneFound``
   0/1, whether or not variant gene is found in HGMD.

``hgmdVarFound``
   0/1, whether or not variant itself is found in HGMD.

``curationScoreHGMD``
   1/2/3 (Low/Medium/High), curated with ``hgmdGeneFound`` and ``hgmdVarFound``.

``omimGeneFound``
   0/1, whether or not variant gene is found in OMIM.

``omimVarFound``
   0/1, whether or not variant itself is found in OMIM.

``curationScoreOMIM``
   1/2/3 (Low/Medium/High), curated with ``omimGeneFound`` and ``omimVarFound``.

``dominant``
   0/1. Whether the variant gene is annotated as dominant in OMIM.

``recessive``
   0/1. Whether the variant gene is annotated as recessive in OMIM

``hgmd_rs``
   HGMD rank score, interpreted as relative probabilities of pathogenicity.

``c_ClinVar_*``
   Expansions of variant annotation from ClinVar. One-hot encoded.

``c_CLNREVSTAT``
   The ClinVar Review status for the same protein change in ClinVar

``c_HGMD_Exp_*``
   Expansions of variant annotation from HGMD. One-hot encoded.

``c_isBLB``
   The original variant is annotated as Benign in ClinVar

``c_isPLP``
   The original variant is annotated as Pathogenic or likely pathogenic in ClinVar

``c_RANKSCORE``
   The HGMD RANKSCORE adapted from the original HGMD database

``nc_ClinVar_Exp``
   Non-coding variant expansion (2bp upstream or downstream of the original variants position)

``nc_CLNREVSTAT``
   non-coding variant expansion (2bp upstream or downstream of the original variants position)

``nc_HGMD_Exp``
   non-coding variant expansion (2bp upstream or downstream of the original variants position)

``nc_isBLB``
   The original variant is annotated as Benign in ClinVar

``nc_isPLP``
   The original variant is annotated as Pathogenic or likely pathogenic in ClinVar

``nc_RANKSCORE``


**Variant Impact**

``cons_*``
   Variant consequence type is one-hot encoded.
   Complete list: 
               'transcript_ablation', 'splice_acceptor_variant', 'splice_donor_variant', 'stop_gained', 
               'frameshift_variant', 'stop_lost', 'start_lost', 'transcript_amplification', 'inframe_insertion', 
               'inframe_deletion', 'missense_variant', 'protein_altering_variant', 'splice_region_variant',
               'splice_donor_5th_base_variant', 'splice_donor_region_variant'

``IMPACT``
   Integer 0-4 (None, Modifier, Low, Moderate, High). Subjective impact classification of consequence type.

``IMPACT.from.Tier``


**In Silico Prediction**

``CADD_phred``
   CADD Phred score

``DANN_score``
   DANN score

``fathmm_MKL_coding_score``
   fathmm-MKL coding socre from dbNSFP

``FATHMM_score``
   FATHMM score from dbNSFP, minimum value selected.

``M_CAP_score``
   M-CAP score

``MutationAssessor_score``
   MutationAssessor score, maximum value selected.

``Polyphen2_HDIV_score``
   Polyphen2 HDIV score, maximum value selected.

``Polyphen2_HVAR_score``
   Polyphen2 HVAR score, maximum value selected.

``REVEL_score``
   REVEL score, maximum value selected.

``SIFT_score``
   SIFT score, minimum value selected.

**Inferred Inheritance**

``No.Var.H``
   Gene level, Number of High IMPACT variants in the patient for candidate gene

``No.Var.HM``
   Gene level, Number of High or Moderate IMPACT variants in the patient for candidate gene

``No.Var.L``
   Gene level, Number of Low IMPACT variants in the patient for candidate gene

``No.Var.M``
   Gene level, Number of Moderate IMPACT variants in the patient for candidate gene

``TierAD``
   1~4, Dominant Inheritance Score. The lower the more pathogenic

``TierAR``
   1~4, Recessive Inheritance Score. The lower the more pathogenic

``TierAR.adj``
   1~4, Adjusted Recessive Inheritance Score. For a candidate gene, if a rare intronic variant observed together with a high IMPACT variant, adjusted

``AD.matched``
   0/1, ``TierAD <= 2`` and  ``dominant == 1``

``AR.matched``
   0/1, ``TierAR <= 2`` and ``recessive == 1``

``zyg``
   Variant zygosity, 1: heterozygous, 2: homozygous.

**Minor Allele Frequency**

``ESP6500_AA_AF``
   ESP6500 African American Allele Frequency

``ESP6500_EA_AF``
   ESP6500 European American Allele Frequency

``gnomadAF``
   gnomAD exome Allele Frequency

``gnomadAFg``
   gnomAD genome Allele Frequency

**Phenotype Matching**

``clinVarSymMatchFlag``
   0/1, whether OMIM variant phenotype matches condition in ClinVar.

``hgmdSymptomSimScore``
   Similarity score between patient phenotype and variant phenotype in HGMD.

``hgmdSymMatchFlag``
   0/1, whether ``hgmdSymptomSimScore >= 0.2``

``omimSymptomSimScore``
   Similarity score between patient phenotype and variant phenotype in OMIM.

``omimSymMatchFlag``
   0/1, whether ``omimSymptomSimScore >= 0.2``

``phrank``
   Phrank measures phenotype sets similarity of the patient phenotype with phenotype linked to a candidate gene.

``diffuse_Phrank_STRING``
   A phenotype score is derived through network diffusion, utilizing the String network and employing the Phrank score as the initial seed score.

**others**

``simple_repeat``
   0/1, whether variant is in simple repeat regions.

``spliceAImax``
   Maximum of SpliceAI score among DS_AG, DS_AL, DS_DG, and DS_DL.


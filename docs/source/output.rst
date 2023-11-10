.. _output:

*********************
AIM Output Details
*********************


Outout Files
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

Outout Explanation
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


Variant Annotation / AIM Features
-------------------------------------

AIM takes information from various database and perform feature engineering, to make the final prediction.
Here, we provide both raw and engineered features as part of the output.

In making the decision, AIM takes various types of information into consideration.

**Conservation**

``GERPpp_RS``
   GERP++_RS, a score indicating the level of evolutionary constraint at a specific genomic position, with positive values suggesting functional conservation and negative values neutrality.

``LRT_Omega``

``LRT_score``

``phyloP100way_vertebrate``
   Conservation score from the phyloP method, with positive scores indicating conservation and negative scores indicating fast evolution.

**Constrain**

``hom``

``decipherVarFound``
   0/1. Whether the variant is found in DECIPHER database. 

``dgvVarFound``
   0/1. Whether the variant is found in DGV database.

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

``c_ClinVar_*``

``c_CLNREVSTAT``

``c_HGMD_Exp_*``

``c_isBLB``

``c_isPLP``

``c_RANKSCORE``

``CLASS``

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

``recessive``

``hgmd_rs``
   HGMD rank score, interpreted as relative probabilities of pathogenicity.

``nc_ClinVar_Exp``

``nc_CLNREVSTAT``

``nc_HGMD_Exp``

``nc_isBLB``

``nc_isPLP``

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

``AD.matched``

``AR.matched``

``No.Var.H``

``No.Var.HM``

``No.Var.L``

``No.Var.M``

``TierAD``

``TierAR``

``TierAR.adj``

``zyg``

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

``diffuse_Phrank_STRING``

**others**

``simple_repeat``
   0/1, whether variant is in simple repeat regions.

``spliceAImax``
   Maximum of SpliceAI score among DS_AG, DS_AL, DS_DG, and DS_DL.


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

``LRT_Omega``

``LRT_score``

``phyloP100way_vertebrate``

**Constrain**

``conservationScoreDGV``

``conservationScoreGnomad``

``conservationScoreOELof``

``decipherVarFound``

``dgvVarFound``

``gnomadGeneOELof``

``gnomadGeneOELofUpper``

``gnomadGenePLI``

``gnomadGeneZscore``

``hom``

**Disease Database**

``c_ClinVar_*``

``c_CLNREVSTAT``

``c_HGMD_Exp_*``

``c_isBLB``

``c_isPLP``

``c_RANKSCORE``

``CLASS``

``clinVarGeneFound``

``clinvarNumB``

``clinvarNumLB``

``clinvarNumLP``

``clinvarNumP``

``clinVarVarFound``

``curationScoreClinVar``

``curationScoreHGMD``

``curationScoreOMIM``

``dominant``

``hgmd_rs``

``hgmdGeneFound``

``hgmdVarFound``

``isB/LB``

``isP/LP``

``nc_ClinVar_Exp``

``nc_CLNREVSTAT``

``nc_HGMD_Exp``

``nc_isBLB``

``nc_isPLP``

``nc_RANKSCORE``

``omimGeneFound``

``omimVarFound``

``recessive``


**Variant Impact**

``cons_*``

``IMPACT``

``IMPACT.from.Tier``


**In Silico Prediction**

``CADD_phred``

``DANN_score``

``fathmm_MKL_coding_score``

``FATHMM_score``

``M_CAP_score``

``MutationAssessor_score``

``Polyphen2_HDIV_score``

``Polyphen2_HVAR_score``

``REVEL_score``

``SIFT_score``

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

``ESP6500_EA_AF``

``gnomadAF``

``gnomadAFg``


**Phenotype Matching**

``clinVarSymMatchFlag``

``hgmdSymMatchFlag``

``hgmdSymptomScore``

``hgmdSymptomSimScore``

``omimSymMatchFlag``

``omimSymptomSimScore``

``phrank``

``diffuse_Phrank_STRING``

**others**

``simple_repeat``

``spliceAImax``


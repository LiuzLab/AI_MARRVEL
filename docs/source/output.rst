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

In making the decision, AIM takes 6 modules of information into consideration.

**Disease database matching**

**Variant Conservation**

**Inferred Inheritance pattern**

**Gene network diffusion**


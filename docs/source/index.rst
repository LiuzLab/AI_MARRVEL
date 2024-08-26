******************************************
Welcome to AI-MARRVEL documentation!
******************************************


**AI-MARRVEL (AIM)** is an AI system for rare genetic disease diagnosis.  

It takes as input the patient VCF file and the phenotype (formatted in HPO) to predict the causal variant(s).   
In making the prediction, it takes variant annotation from `MARRVEL <https://marrvel.org/>`_ database and more, 
and generates **prediction score** + **confidence score** as output.

Web Interface
=====================

You can use AI-MARRVEL from our `website <https://ai.marrvel.org/>`_ or follow the documentation to run locally.


Prerequisites
=====================

* Unix-like operating system
* Docker installed
* At least 480G space to download the preprocessed database
* 16G free memeory recommended for whole-exome sequencing data, greater memory is required for whole-genome sequencing data.



Citation
=====================

Mao, D., Liu, C., Wang, L., AI-Ouran, R., Deisseroth, C., Pasupuleti, S., Kim, S. Y., Li, L., Rosenfeld, J. A., Meng, L., B
urrage, L. C., Wangler, M. F., Yamamoto, S., Santana, M., Perez, V., Shukla, P., Eng, C. M., Lee, B., Yuan, B., … Liu, Z. (2024). 
**AI-MARRVEL—A knowledge-driven ai system for diagnosing mendelian disorders. NEJM AI, 1(5).** <a href="https://doi.org/10.1056/AIoa2300009">https://doi.org/10.1056/AIoa2300009</a>


License
=====================
AI-MARRVEL is licensed under `GPL-3.0 <https://github.com/LiuzLab/AI_MARRVEL/blob/main/LICENSE>`_.

© Copyright 2024 by Zhandong Liu's lab at Baylor College of Medicine.

.. toctree::
   :hidden:

   installation
   runningaim
   input
   output
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



.. note::

   AI-MARRVEL and web interface are under active development.



Prerequisites
=====================

* Unix-like operating system
* Docker installed
* At least 480G space to download the preprocessed database
* 32G free memeory recommended, but it can be smaller by passing memory limitation argument



Citation
=====================
**Enhancing Molecular Diagnostics of Mendelian Disorders with AI-MARRVEL: A Knowledge-Driven Artificial Intelligence Approach.**
Click `here <https://dx.doi.org/10.2139/ssrn.4465963>`_ for details.


License
=====================
AI-MARRVEL is licensed under `GPL-3.0 <https://github.com/LiuzLab/AI_MARRVEL/blob/main/LICENSE>`_.

© Copyright 2023 by Zhandong Liu's lab at Baylor College of Medicine.

.. toctree::
   :hidden:

   installation
   runningaim
   input
   output
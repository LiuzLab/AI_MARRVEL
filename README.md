<p align="center">
  <img width="50%" align="center" src="https://raw.githubusercontent.com/LiuzLab/AI_MARRVEL/main/docs/images/logo.v1.png" alt="logo">
</p>
  <h1 align="center">
  AI-MARRVEL
</h1>

<p align="center">
  <a style="text-decoration:none">
    <img src="https://img.shields.io/badge/AI_MARRVEL-v1.0.1-blue.svg"/>
  </a>
  <a href='https://ai-marrvel.readthedocs.io/en/latest/?badge=latest'>
      <img src='https://readthedocs.org/projects/ai-marrvel/badge/?version=latest' alt='Documentation Status' />
  </a>
  <a style="text-decoration:none">
    <img src="http://img.shields.io/badge/Preprint-DOI:10.2139/ssrn.4465963-lightgreen.svg"/>
  </a>
</p>

**AI-MARRVEL (AIM)** is an AI system for rare genetic disease diagnosis.  

It takes as input patient VCF and phenotype (formatted with HPO) to predict the causal variant(s).    
In making prediction, it takes variant annotation from [MARRVEL](https://marrvel.org/) database and more, 
and generates **prediction score** + **confidence score** as output.  


You can use AI-MARRVEL from our [website](https://ai.marrvel.org/) or follow the [documentation](https://ai-marrvel.readthedocs.io/en/latest/) to run locally.



## Quick Start

### Install Required Data Dependencies
AIM utilizes various databases for variant annotation, all of which have been compiled and are available for download. We use AWS S3 for data access, and the data can be downloaded by following these steps:

1. **Install the AWS CLI**: Follow the instructions provided in the [AWS CLI Installation Guide](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html).
2. **Navigate to Your Desired Directory**: Change to the directory where you want your data dependencies downloaded. For example, in Ubuntu, use:
   ```bash
   cd <desired/folder/path>
3. Use the following command to sync the S3 bucket to your local directory:  
    ``` bash
    aws s3 sync s3://aim-data-dependencies-public . --no-sign-request
### Get the software
AIM is released as a Docker image for easy distribution. To get it:
```
docker pull chaozhongliu/aim-lite:latest
```

### Run with your sample
```
nextflow run .  --ref_dir <PATH_TO_REFERENCE_DIRECTORY>
                --input_vcf <PATH_TO_INPUT_VCF_FILE>
                --input_hpo <PATH_TO_INPUT_HPO_FILE>
                --outdir <PATH_TO_OUTPUT_DIRECTORY>
                --run_id [Sample Id] # Optional, default: 1
                --ref_ver [Reference genome: hg19/hg38] # Optional, default: hg19
```

Alternatively, the pipeline can be executed with a parameter file (yaml)
```
nextflow run . -params-file params.yaml
```
NOTE: You need to create `params.yaml` by copying [params.yaml.example](params.yaml.example) file and follow the instruction.

## License
AI-MARRVEL is licensed under GPL-3.0. You are welcomed to use it for research purpose.  
For business purpose, please contact us for licensing.


## Disclaimer
- Some of the data and software included in the distribution may be subject to third-party constraints. Users of the data and software are solely responsible for establishing the nature of and complying with any such restrictions.
- AI-MARRVEL provides this data and software in good faith, but make no warranty, express or implied, nor assume any legal liability or responsibility for any purpose for which they are used.

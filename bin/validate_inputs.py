#!/usr/bin/env python3.8
import sys
import os
import json
import csv

def validate_file_type(file_path, extensions, error_msg_no_exist, error_msg_wrong_type):
    # Check if the path exists and is a file
    if not os.path.isfile(file_path):
        return False, error_msg_no_exist
    
    # Check if the file has the correct extension
    correct_extension = any(file_path.endswith(ext) for ext in extensions)
    # print(f"Resolved VCF file path: {file_path}, correct_extension: {correct_extension}")
    
    if not correct_extension:
        return False, error_msg_wrong_type
    
    return True, f"File exists and is of the correct type: {file_path}"

def validate(user_inputs, vcf, hpo, refDir):
    errors = []
    missing_params = []
    results = []
    required_user_inputs = [ 'input_vcf', 'input_hpo', 'ref_dir', 'run_id', 'ref_ver']

    # Check Required User Inputs for missing parameters
    for param in required_user_inputs:
        if param not in user_inputs or user_inputs[param] is None or str(user_inputs[param]).strip() == "":
            missing_params.append(param)

    if len(missing_params) > 0:
        error_message = (
            "The following required parameters are missing: "
            f"{', '.join(missing_params)} \n"
            "Please provide values for these parameters and run the pipeline again."
        )
        results.append({
            "parameter": "Missing Parameters",
            "status": "Invalid",
            "error_message": error_message,
            "error_type": "Missing_Parameters_Error"
        })
        return results
    else:
        # print(f"Resolved VCF file path: {os.path.isfile(vcf)}")
        # print(f"Resolved HPO file path: {os.path.isfile(hpo)}")
        # print(f"Resolved Reference Directory path: {os.path.isdir(refDir)}")

        # Checking Required User Parameters
        # VCF
        vcf_error_no_exist = f"""
            Input VCF file '{user_inputs.get("input_vcf")}' provided in --input_vcf parameter does not exist.\n 
            Please make sure VCF file exists in the path and has the right format. 
        """
        vcf_error_wrong_type = f"""
            Input VCF file '{user_inputs.get("input_vcf")}' provided in --input_vcf parameter is not of the correct type.\n 
            Please provide a valid VCF file with either ".vcf" or ".vcf.gz" extension.
        """
        is_valid_vcf, vcf_message = validate_file_type(vcf, [".vcf", ".vcf.gz"], vcf_error_no_exist, vcf_error_wrong_type)
        if not is_valid_vcf:
            errors.append(vcf_message)

        # HPO
        hpo_error_no_exist = f"""
            Input HPO file '{user_inputs.get("input_hpo")}' provided in --input_hpo parameter does not exist.\n
            Please make sure HPO file exists in the path and has the right format. 
        """
        hpo_error_wrong_type = f"""
            Input HPO file '{user_inputs.get("input_hpo")}' provided in --input_hpo parameter is not of the correct type.\n
            Please provide a valid HPO file with the ".txt" extension.
        """
        is_valid_hpo, hpo_message = validate_file_type(hpo, [".txt"], hpo_error_no_exist, hpo_error_wrong_type)
        if not is_valid_hpo:
            errors.append(hpo_message)


        # Data Dependencies      
        if not os.path.isdir(refDir):
            errors.append(f"""
                          Reference Data dependence directory '{user_inputs.get("ref_dir")}' provided in --ref_dir parameter does not exist.\n
                          Please make sure Data dependence exists in the path and has the right format.
                          """)

        # Ref Genome 
        ref_ver = user_inputs.get("ref_ver", "").strip()
        if ref_ver not in ['hg19', 'hg38']:
            errors.append(f"Invalid reference version '{ref_ver}'. Please enter as either '--ref_ver hg19' or '--ref_ver hg38'.")

        ## TODO : Add remaining validation to this file here. 
        

        # collecting errors messages to results
        if errors:
            for error in errors:
                results.append({
                    "parameter": "Validation Check",
                    "status": "Invalid",
                    "error_message": error,
                    "error_type": "Input_Validation_Error"
                })
        else:
            results.append({
                "parameter": "Validation Check",
                "status": "Valid",
                "error_message": "",
                "error_type": "No_Error"
            })
        return results

if __name__ == "__main__":
    user_inputs = json.loads(sys.argv[1])
    input_vcf = sys.argv[2]
    input_hpo = sys.argv[3]
    input_ref = sys.argv[4]

    # Perform validation
    results = validate(user_inputs,input_vcf,input_hpo,input_ref)

    # Write the results to a CSV file
    with open('validation_results.csv', 'w', newline='') as csvfile:
        fieldnames = ['parameter', 'status', 'error_message', 'error_type']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        for result in results:
            writer.writerow(result)

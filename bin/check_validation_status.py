#!/usr/bin/env python3.8
import sys
import csv
import os
import textwrap

def check_validation(file_path):
    # Verify if the file is accessible
    if not os.path.isfile(file_path):
        print(f"ERROR: File not found: {file_path}")
        sys.exit(1)
    else:
        print("Validation file found")

    # Read the CSV file
    invalid_status = []
    try:
        with open(file_path, 'r') as csv_file:
            reader = csv.DictReader(csv_file)
            
            # Trim whitespace from headers
            reader.fieldnames = [header.strip() for header in reader.fieldnames]
            
            invalid_status = []
            for row in reader:
                status_value = row.get('status', '').strip().lower()
                if status_value == "invalid":
                    invalid_status.append(row)
        
        if invalid_status:
            print("ERROR: Invalid inputs:")
            for row in invalid_status:
                print(textwrap.dedent(f"""\
                    -------------------------------------------------------
                    Type    : {row.get('error_type', '').strip()}
                    Message : {row.get('error_message', '').strip()}
                    -------------------------------------------------------
                """))
            print("To see usage and available parameters run `nextflow run main.nf --help`")
            sys.exit(128)
        
        return True

    except KeyError as e:
        print(f"ERROR: KeyError - the key '{e.args[0]}' was not found in the CSV file.")
        sys.exit(1)

    except Exception as e:
        print(f"ERROR: extracting from csv file: {file_path}, Exception: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    csv_file_path = sys.argv[1]
    validation_status = check_validation(csv_file_path)
    print(validation_status)

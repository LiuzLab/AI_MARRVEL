##PEDIGREE=<Child=1000021,Mother=1000023,Father=1000022>
import gzip
from glob import glob
import sys
import argparse



def get_PED(file_name):
    ped_line = ''
    with gzip.open(file_name, 'rb') as F:
        for line in F:
            if not line.startswith(b'##'):
                break
            if line.startswith(b'##PEDIGREE'):
                ped_line = line
                break
            else:
                continue
    if ped_line == '':
        return None

    ped_line = ped_line.decode("utf-8")
    
    line = ped_line.lstrip('##PEDIGREE=<')
    line = line.rstrip('>\n')
    if 'Child' in line and 'Mother' in line and 'Father' in line:
        line = line.split(',')
        line = [(i.split('=')[0], i.split('=')[1]) for i in line]
        ped = {i[0]:i[1] for i in line}
        return ped
    else:
        return None


def make_PED(ped,out_dir):
    # Example: 
    # 1 2020717 2020721 2020718 0 2
    # 1 2020718 0 0 2 1
    # 1 2020721 0 0 1 1
    ped_file = open(f"{out_dir}/PEDs/{ped['Child']}.ped", 'w')
    ped_file.write(f"""1\t{ped['Child']}\t{ped['Father']}\t{ped['Mother']}\t0\t2\n1\t{ped['Mother']}\t0\t0\t2\t1\n1\t{ped['Father']}\t0\t0\t1\t1""")
    ped_file.close()



def generate_new_cohort(vcf_path, hpo_path, out_dir, out_prefix=''):

    file_names = glob(f"{vcf_path}/*")
    file_names = [i for i in file_names if 'gz' in i]
    hpo_files = glob(hpo_path)
    if out_prefix == '':
        out_prefix = 'sample'
    samples_list = open(f"{out_dir}/{out_prefix}.txt", 'w')

    #complete_samples = []
    no_ped_samples = []
    no_hpo_samples = []

    for file_name in file_names:
        ped = get_PED(file_name)

        if ped is not None:
            if not glob(f"{hpo_path}/*{ped['Child']}*"):
                no_hpo_samples.append(file_name)
                continue

            hpo_file = glob(f"{hpo_path}/*{ped['Child']}*")[0]
            make_PED(ped, out_dir)
            samples_list.write(f"{ped['Child']}\t{ped['Mother']}\t{ped['Father']}\t{file_name}\t{hpo_file}\n")
            #complete_samples.append(file_name)

        else:
            no_ped_samples.append(file_name)
            continue

    samples_list.close()

    miss_ped = open(f'{out_dir}/missing_ped_{out_prefix}.txt', 'w')
    miss_ped.write('\n'.join(no_ped_samples))
    miss_ped.close()

    miss_hpo = open(f'{out_dir}/missing_hpo_{out_prefix}.txt', 'w')
    miss_hpo.write('\n'.join(no_hpo_samples))
    miss_hpo.close()



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-varPath", "--varPath",  help="Joint VCF file folder full path")
    parser.add_argument("-hpoPath", "--hpoPath",  help="HPO file folder full path")
    parser.add_argument("-outDir", "--outDir",  help="Output direction full path")
    parser.add_argument("-prefix", "--prefix",  default='', help="Output files prefix")
    
    args = parser.parse_args()
    generate_new_cohort(args.varPath, args.hpoPath, args.outDir, args.prefix)



main()










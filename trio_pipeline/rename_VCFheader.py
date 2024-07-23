import sys

proband_id = sys.argv[1]
paternal_id = sys.argv[2]
maternal_id = sys.argv[3]
header_file = sys.argv[4]

sample_list = [proband_id, paternal_id, maternal_id]

header = open(header_file, 'r').readlines()[0]
header = header.strip('\n').split('\t')[-3:]
newheader = []

for sample in header:
	newheader += [i for i in sample_list if i in sample]

assert len(newheader) == 3

sample_name = open("/out/sample_ids.txt", 'w')
sample_name.write('\n'.join(newheader))
sample_name.close()





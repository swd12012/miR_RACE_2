import os
import shutil
import sys
import re
import numpy as np
import pandas as pd

# Find all the .gz files and write the output to ATAClist after splitting output into lines

dir1 = '/dfs6/pub/swdu/microRNA/20210610_miR_RACE/data/raw/'
os.chdir(dir1)

output = os.popen('find *.gz').read()
RACE_list = output.split()
print(RACE_list)

#Remove unnecessary file names by string substitution


sample_names = []

for i in RACE_list:
	i = re.sub('(mR191-L1-)', '', i)
	i = re.sub('-Sequences.txt.gz', '.fq.gz', i)
	i = re.sub('[ATCG]{8}-[ATCG]{8}-', '', i)
	sample_names.append(i)

#print(sample_names)

#Create dictionary linking old names to new names
RACE_dic = {RACE_list[i]: sample_names[i] for i in range(len(RACE_list))}
#print(RACE_dic)

for key in RACE_dic:
	os.rename(key, RACE_dic[key])

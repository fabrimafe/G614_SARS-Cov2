import argparse as ap
from random import shuffle
from datetime import datetime

parser = ap.ArgumentParser()

parser.add_argument('-i', '--infile', help='Provide a nexus file (but should work also with a phylip format)', required=True, type=str)
parser.add_argument('-d', '--last_date', help='Date of the most recent sequence as YYYY-MM-DD', required=True, type=str)
parser.add_argument('-f', '--minfreq', help='Minimum frequency of minor allele', required=True, type=float)
args = parser.parse_args()

#launch the function
input_file = args.infile 
last_date = args.last_date
last = datetime.strptime(last_date, "%Y-%m-%d")
minfreq = args.minfreq

output_file = input_file.replace('.nex','')+'.varfreqmin'+str(minfreq)+'.txt'

seq_list = []
sample_list = []
date_list = []
countries_list = []
subcountries_list = []

with open(input_file, 'r') as handle:
	for line in handle:
		if len(line) < 1000:
			continue
		elif 'env' in line:
			continue
		elif 'reference' in line:
			continue
		elif 'hCoV' in line:
			row = line.rstrip().split('|')
			country = row[0].split('/')[1].replace(' ','')
			subcountry = row[0].split('/')[2][0:2]
			ID = row[1]
			day = row[2].split()[0]
			if len(day.split('-')) == 3:
				date_seq = datetime.strptime(day, "%Y-%m-%d")
				delta = last - date_seq
				days_BP = delta.days
				seq = row[2].split()[1].replace('n','?').replace('-','?')
				seq_list.append(seq)
				sample_list.append(ID)
				date_list.append(days_BP)
				countries_list.append(country)
				subcountries_list.append(subcountry)


tot_seq = len(seq_list)

nuc_list = list(zip(*seq_list))

variant_list = []
position_list = []
for j in range(0,len(seq_list[0])):
	data = [i for i in nuc_list[j] if i != '?']
	alleles = set(data)
	count_alleles = len(alleles)
	if len(data) < tot_seq*0.75:
		continue
	elif count_alleles < 2:
		continue
	else:
		freqs = []
		for i in list(alleles):
			freqs.append(data.count(i))
		print(alleles, freqs, sorted(freqs)[-1]/len(data))
		if sorted(freqs)[-1]/len(data) < 1-minfreq:
			variant_list.append(nuc_list[j])
			position_list.append(j+1)


variants = list(zip(*variant_list))

print('Variants: '+str(len(position_list)))

china_pref_names = ['Beijing',
'Chongqing',
'Changzhou',
'Foshan',
'Fujian',
'Fuyang',
'Fuzhou',
'Ganzhou',
'Guangdong',
'Guangzhou',
'Harbin',
'Hangzhou',
'Hefei',
'Jian',
'Jiangsu',
'Jiangxi',
'Jingzhou',
'Jiujiang',
'Liaoning',
'Lishui',
'Nanchang',
'NanChang',
'Pingxiang',
'Shandong',
'Shanghai',
'Shangrao',
'Shaoxing',
'Shenzhen',
'Sichuan',
'Tianmen',
'Wuhan',
'Wuhan-Hu-1',
'Xinyu',
'Yichun',
'Yingtan',
'Yunnan',
'Zhejiang']


out = open(output_file, 'w')

out.write('ID\tdaysbefore'+last_date+'\tcountry\tsubcountry')
for i in position_list:
	out.write('\t')
	out.write(str(i))
out.write('\n')


for i in range(0,len(seq_list)):
	if countries_list[i] in china_pref_names:
		out.write(sample_list[i]+'\t'+str(date_list[i])+'\tChina\t'+countries_list[i]+'\t'+'\t'.join(variants[i])+'\n')
	else:
		out.write(sample_list[i]+'\t'+str(date_list[i])+'\t'+countries_list[i]+'\t'+subcountries_list[i]+'\t'+'\t'.join(variants[i])+'\n')
out.close()





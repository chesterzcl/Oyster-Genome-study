import numpy as np
import pandas as pd


class intvlist:
	def __init__(self):
		self.intervals = []
    
	def insert(self, new_interval):
		start, end = new_interval
		merged = []
		inserted = False
        
		for interval in self.intervals:
			if interval[1] < start:
				merged.append(interval)
			elif end < interval[0]:
				if not inserted:
					merged.append((start, end))
					inserted = True
				merged.append(interval)
			else:
				start = min(start, interval[0])
				end = max(end, interval[1])
                
		if not inserted:
			merged.append((start, end))
        
		self.intervals = sorted(merged, key=lambda x: x[0])
    
	def total_length(self):
		return sum(end - start for start, end in self.intervals)
	
	def get_intervals(self):
		return self.intervals


dir="/Users/lizicheng/Desktop/Data/dog_genome/oyster"
file="CV_primary_polished_scfd.coords"

line_cter=0

chr_dict={}
scfd_dict={}
chr_total_len_dict={}
scfd_total_len_dict={}


with open(dir+"/"+file,'r') as file:
	for line in file:
		line_cter+=1
		line_list=line.strip().split('|')
		if line_cter>5:
			# len_info=line_list[2]
			# len_info_chr=len_info.strip().split()[0]
			# len_info_scfd=len_info.strip().split()[1]
			chr_intv=line_list[0]
			chr_intv_start=int(chr_intv.strip().split()[0])
			chr_intv_end=int(chr_intv.strip().split()[1])

			total_len_info=line_list[4]
			total_chr_len=total_len_info.strip().split()[0]
			total_scfd_len=total_len_info.strip().split()[1]

			scfd_info=line_list[-1]
			chr=scfd_info.strip().split()[0]
			scfd=scfd_info.strip().split()[1]

			if chr not in chr_dict:
				chr_dict[chr]={}
			else:
				if scfd not in chr_dict[chr]:
					chr_dict[chr][scfd]=intvlist()
				chr_dict[chr][scfd].insert((chr_intv_start,chr_intv_end))

			# if scfd not in scfd_dict:
			# 	scfd_dict[scfd]={}
			# else:
			# 	if chr not in scfd_dict[scfd]:
			# 		scfd_dict[scfd][chr]=0
			# 	scfd_dict[scfd][chr]+=int(len_info_scfd)

			if chr not in chr_total_len_dict:
				chr_total_len_dict[chr]=int(total_chr_len)
			if scfd not in scfd_total_len_dict:
				scfd_total_len_dict[scfd]=int(total_scfd_len)


chr_dict_frac={}
chr_dict_frac_sorted={}

del chr_dict['NC_007175.2']

for key_chr in chr_dict:
	chr_dict_frac[key_chr]={}
	for key_scfd in chr_dict[key_chr]:
		chr_dict_frac[key_chr][key_scfd]=chr_dict[key_chr][key_scfd].total_length()/chr_total_len_dict[key_chr]
	chr_dict_frac_sorted[key_chr]=dict(sorted(chr_dict_frac[key_chr].items(),key=lambda item: item[1],reverse=True))


print(chr_total_len_dict)
print(chr_dict_frac_sorted)

for key_chr in chr_dict_frac_sorted:
	print(key_chr,chr_total_len_dict[key_chr])
	cter=0
	for key_scfd in chr_dict_frac_sorted[key_chr]:
		if cter<3:
			print(key_scfd,chr_dict_frac_sorted[key_chr][key_scfd])
		cter+=1

		
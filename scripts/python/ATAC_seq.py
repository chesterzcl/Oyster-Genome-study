import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import math
import pandas as pd


chr_len_dict={
"NC_035780.1":65668440,
"NC_035781.1":61752955,
"NC_035782.1":77061148,
"NC_035783.1":59691872,
"NC_035784.1":98698416,
"NC_035785.1":51258098,
"NC_035786.1":57830854,
"NC_035787.1":75944018,
"NC_035788.1":104168038,
"NC_035789.1":32650045
}

mapped_reads_dict={
"3L":22815331,
"4L":21409280,
"5L":18987022,
"6L":23646885,
"7L":23356852,
"8L":21753757,
"9L":1134953,
"10L":1369823,
"11L":19907575,
"12L":20480122
}

dir="/Users/lizicheng/Desktop/Data/dog_genome/oyster/ATAC_L_peak"
for i in range(3,13):
	sample_name=str(i)+"L"
	file_name=sample_name+"_ATAC_st_filtered_shifted_mtrm_st_peak_peaks.narrowPeak"

	peak_df=pd.read_csv(dir+"/"+file_name,sep='\t',header=None)
	scaling_fact=20
	total_reads=mapped_reads_dict[sample_name]
	insert_size_list_allchr=[]

	for chr in chr_len_dict:
		peak_list_nfr=[[],[],[]]
		peak_list_mono=[[],[],[]]
		peak_list_di=[[],[],[]]
		peak_list_multi=[[],[],[]]

		#Peak distribution across chromosome
		for i,j in peak_df.iterrows():
			if j[0]==chr:
				peak_size=j[2]-j[1]
				peak_pos=j[1]+j[9]
				peak_str=j[6]/scaling_fact
				if peak_size<=100:
					peak_list_nfr[0].append(peak_pos)
					peak_list_nfr[1].append(peak_str)
					peak_list_nfr[2].append(peak_size)
				elif peak_size>=120 and peak_size<=180:
					peak_list_mono[0].append(peak_pos)
					peak_list_mono[1].append(peak_str)
					peak_list_mono[2].append(peak_size)
				elif peak_size>=250 and peak_size<=350:
					peak_list_di[0].append(peak_pos)
					peak_list_di[1].append(peak_str)
					peak_list_di[2].append(peak_size)

		fig,ax=plt.subplots(figsize=(10,4))

		peak_list_dict={"Nucleosome-free":peak_list_nfr,"Mono-nucleosomal":peak_list_mono,"Di-nucleosomal":peak_list_di}

		for key,value in peak_list_dict.items():

			for i in range(len(value[0])):
				rect=patches.Rectangle((value[0][i],0),value[2][i],value[1][i],color='black',alpha=0.7)
				ax.add_patch(rect)

			ax.set_xlim(1,chr_len_dict[chr])
			ax.set_ylim(0,10)
			ax.set_xlabel(f"Genomic Position on Chr {chr} ")
			ax.set_ylabel("Peaks")
			ax.set_title(f"{key} Peaks in Chr {chr} for sample {sample_name}")
			plt.savefig(dir+"/"+sample_name+"_"+key+"_"+chr+".png") 

	
	#Length distribution of accessible regions
	# plt.close('all')
	# plt.figure(figsize=(8, 6)) 
	# plt.hist(insert_size_list_allchr,bins=10,color='skyblue', edgecolor='black')
	# plt.xlabel("Peak Width (bp)")
	# plt.ylabel("Number of Peaks")
	# plt.title("Distribution of ATAC‑seq Peak Width for sample "+sample_name)
	# plt.show()

		#Strength distribution of peaks

		# plt.hist(peak_str_lst,bins=50,color='skyblue', edgecolor='black')
		# plt.xlabel("Peak Signal Value")
		# plt.ylabel("Frequency")
		# plt.title("Distribution of ATAC‑seq Peak Signal Values")
		# plt.show()














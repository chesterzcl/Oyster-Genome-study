import numpy as np
from scipy.stats import fisher_exact

dir="/Users/lizicheng/Desktop/Data/dog_genome/oyster"
ip_file="20cv_df3_hoh_biSNP_filtered_maf05_ann_sig5_sam8_dist.txt"
op_file="20cv_df3_hoh_biSNP_filtered_maf05_ann_sig5_sam8_dist_reformat.txt"

with open(dir+'/'+ip_file, "r") as file,open(dir+"/"+op_file,"w") as o_file:
	for line in file:
		if line[0]!='#':
			line_lst=line.strip().split('\t')
			ref_allele=line_lst[2]
			alt_allele=line_lst[3]
			ref_allele_s=0
			ref_allele_l=0
			alt_allele_s=0
			alt_allele_l=0
			gt_line=""
			dp_line=""
			for i in range(5,len(line_lst)):
				var_list=line_lst[i].split(';')
				gt=var_list[0]
				dp=var_list[1]
				dp_dist=dp.split(',')
				if gt[0]=='0' and gt[-1]=='0':
					if i<=14:
						ref_allele_l=ref_allele_l+2
					else:
						ref_allele_s=ref_allele_s+2
				elif gt[0]=='0' and gt[-1]=='1':
					if i<=14:
						ref_allele_l=ref_allele_l+1
						alt_allele_l=alt_allele_l+1
					else:
						ref_allele_s=ref_allele_s+1
						alt_allele_s=alt_allele_s+1
				elif gt[0]=='1' and gt[-1]=='1':
					if i<=14:
						alt_allele_l=alt_allele_l+2
					else:
						alt_allele_s=alt_allele_s+2


				if gt=='./.':
					dp_line=dp_line+'\t'+'.'+'\t'+'.'
				else:
					dp_line=dp_line+'\t'+dp_dist[0]+'\t'+dp_dist[1]

				gt_line=gt_line+'\t'+ref_allele+'\t'+alt_allele
			oddsr,p=fisher_exact([[alt_allele_l,alt_allele_s],[ref_allele_l,ref_allele_s]], alternative='two-sided')
			logp=-np.log10(p)
			gt_line='\t'+str(ref_allele_l)+'\t'+str(alt_allele_l)+'\t'+str(ref_allele_s)+'\t'+str(alt_allele_s)+'\t'+str(logp)+gt_line
			dp_line='\t'+'\t'+'\t'+'\t'+'\t'+dp_line
			gt_line=line_lst[0]+'\t'+line_lst[1]+'\t'+line_lst[2]+'\t'+line_lst[3]+'\t'+line_lst[4]+gt_line
			dp_line='\t'+'\t'+'\t'+'\t'+dp_line
			
			o_file.write(gt_line+'\n')
			o_file.write(dp_line+'\n')
		else:
			line_lst=line.strip().split('\t')
			header_line=""
			header_line=header_line+line_lst[0][1:]
			for i in range(1,5):
				header_line=header_line+'\t'+line_lst[i]
			header_line=header_line+'\t'+"ref allele in large group"+'\t'+"alt allele in large group"+'\t'+"ref allele in small group"+'\t'+"alt allele in small group"+'\t'+"-log(P_fisher)"
			for i in range(5,len(line_lst)):
				header_line=header_line+'\t'+line_lst[i]+'\t'
			o_file.write(header_line+'\n')



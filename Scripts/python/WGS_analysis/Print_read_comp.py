import pandas as pd
import numpy as np


dir="/Users/lizicheng/Desktop/Data/dog_genome/oyster/extracted_gvcf"

ref_file_name="chr1_sample.fa"

start=10001
end=10100

seq_lst=[]
pos_lst=list(range(start,end+1))
var_lst=[""]*(end-start+1)
dp_lst=[-1]*(end-start+1)
with open(dir+"/"+ref_file_name, 'r') as file:
    for line in file:
    	seq=line.strip()
    	for char in seq:
    		seq_lst.append(char)


seq_df=pd.DataFrame({'pos':pos_lst,'ref':seq_lst,'var':var_lst,'dp':dp_lst})        
print(seq_df)

sample_st={'1S','3S','8S','12S','14S','3L','5L','7L','8L','9L','11L','12L'}

for sample_name in sample_st:

	file_name=sample_name+"_5000.txt"

	var_df=pd.read_csv(dir+"/"+file_name,sep='\t',header=None,comment="#")
	filtered_df = var_df[(var_df[1] >= 10001) & (var_df[1] <= 10100)]
	filtered_df=filtered_df.drop(filtered_df.columns[0], axis=1)
	filtered_df.columns=['pos','name','ref','var','score','info','ann1','ann2','ann3']

	for idx,row in filtered_df.iterrows():
		if row['var']!='<NON_REF>':
			var_all= row['var'].split(',')[0]
			ref=row['ann3'].split(':')[1].split(',')[0]
			alt=row['ann3'].split(':')[1].split(',')[1]
			gt=str(ref)+','+str(alt)
			seq_df.loc[seq_df['pos']==row['pos'],['dp']]=gt
			seq_df.loc[seq_df['pos']==row['pos'],['var']]=row['ref']+","+var_all
			# print(seq_df[seq_df['pos']==row['pos']])

		else:
			depth = row['ann3'].split(':')[1]
			reg_begin=row['pos']
			reg_end=row['ann1'].split('=')[1]
			for i in range(int(reg_begin),int(reg_end)+1):
				if i<=end:
					seq_df.loc[seq_df['pos']==i,['var']]=seq_df[seq_df['pos']==i]['ref']
					seq_df.loc[seq_df['pos']==i,['dp']]=depth
					# print(seq_df[seq_df['pos']==i])

	seq_df.to_csv(dir+'/'+sample_name+"_seq.txt",index=False,sep='\t')





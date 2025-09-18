import numpy as np
from scipy.stats import fisher_exact

dir="/home/zl436/palmer_scratch/vcf/oyster"
vcf="20cv_hoh_biSNP_filtered.vcf"
op_file="20cv_hoh_biSNP_filtered_ann_fe.stats"
op_sig_file="20cv_hoh_biSNP_filtered_sig_dist.txt"
# dir="/Users/lizicheng/Desktop/Data/dog_genome/oyster"
# vcf="20cv_NC_035781.1_ann.vcf"
# op_file="20cv_test.stats"
# op_sig_file="20cv_test_sig.txt"

min_num=8
p_thres=6

# add near the top (after imports)
def is_non_intron_genic(info_field: str) -> bool:
    """
    Return True if ANN contains at least one genic, non-intron effect.
    Accepts exonic/CDS/UTR/noncoding-exon, rejects intron/intergenic/flanks and splice donor/acceptor.
    """
    ann_val = None
    for tag in info_field.split(';'):
        if tag.startswith('ANN='):
            ann_val = tag[4:]
            break
    if not ann_val:  # no annotation → don't keep
        return False

    # Effects we explicitly do NOT want
    disallow = {
        'intron_variant', 'intergenic_region',
        'upstream_gene_variant', 'downstream_gene_variant',
        'splice_acceptor_variant', 'splice_donor_variant'
    }

    # Genic, non-intron effects we DO want
    allow = {
        # coding / protein-changing or neutral
        'missense_variant', 'synonymous_variant', 'stop_gained', 'stop_lost',
        'start_lost', 'start_gained', 'coding_sequence_variant',
        'frameshift_variant', 'inframe_insertion', 'inframe_deletion',
        'protein_altering_variant', 'stop_retained_variant',
        # UTRs
        '5_prime_UTR_variant', '3_prime_UTR_variant',
        # noncoding exons
        'non_coding_transcript_exon_variant',
        # optional: near-exon boundary; exclude if you consider this intronic
        'splice_region_variant'
    }

    # ANN entries are comma-separated; second pipe-delimited field is the effect term
    for entry in ann_val.split(','):
        parts = entry.split('|')
        if len(parts) > 1:
            effect = parts[1]
            if effect in disallow:
                continue
            if effect in allow:
                return True
    return False

header=""
with open(dir+"/"+vcf,'r') as file:
	for line in file:
		if line[0]=='#':
			header=line.strip()

header_ls=header.split('\t')

s_idx_lst=[]
l_idx_lst=[]

for i in range(9,len(header_ls)):
	if header_ls[i][-1]=='L':
		l_idx_lst.append(i)
	else:
		s_idx_lst.append(i)

cter=0

with open(dir+'/'+vcf,'r') as ip,open(dir+'/'+op_file,'w') as op, open(dir+'/'+op_sig_file,'w') as op2:
	
	op.write("#Chromosome\tposition\tref_allele\talt_allele\tannotation\todds_ratio\tp_value\n")
	
	op2_header="#Chromosome\tposition\tref_allele\talt_allele\tannotation"
	
	for idx in l_idx_lst:
		op2_header=op2_header+"\t"+header_ls[idx]
	for idx in s_idx_lst:
		op2_header=op2_header+"\t"+header_ls[idx]
	
	op2_header+="\n"
	print(op2_header)
	op2.write(op2_header)
	# op2.flush()

	for line in ip:
		if line[0]!='#':
			cter+=1
			if cter%100000==0:
				print(str(cter),"variants analyzed.")
			gt_line=line.strip()
			line_lst=gt_line.split('\t')
			if not is_non_intron_genic(line_lst[7]):  # line_lst[7] is the INFO field
				continue
			s_var_count=0
			s_ref_count=0
			s_valid_num=0
			l_var_count=0
			l_ref_count=0
			l_valid_num=0

			for idx in l_idx_lst:
				gt=line_lst[idx].split(':')[0]
				if gt[0]=='0' and gt[-1]=='0':
					l_ref_count+=2
					l_valid_num+=1
				elif gt[0]=='0' and gt[-1]=='1':
					l_var_count+=1
					l_ref_count+=1
					l_valid_num+=1
				elif gt[0]=='1' and gt[-1]=='1':
					l_var_count+=2
					l_valid_num+=1

			for idx in s_idx_lst:
				gt=line_lst[idx].split(':')[0]
				if gt[0]=='0' and gt[-1]=='0':
					s_ref_count+=2
					s_valid_num+=1
				elif gt[0]=='0' and gt[-1]=='1':
					s_var_count+=1
					s_ref_count+=1
					s_valid_num+=1
				elif gt[0]=='1' and gt[-1]=='1':
					s_var_count+=2
					s_valid_num+=1

			if l_valid_num>=min_num and s_valid_num>=min_num:
				oddsr,p=fisher_exact([[l_var_count,s_var_count],[l_ref_count,s_ref_count]], alternative='two-sided')
				ann_str=line_lst[7].split(';')[-1]
				gene_str=ann_str.split('|')[3]
				op_line=line_lst[0]+'\t'+line_lst[1]+'\t'+line_lst[3]+'\t'+line_lst[4]+'\t'+gene_str+'\t'+str(oddsr)+'\t'+str(p)+'\n'
				op.write(op_line)

				logp=-np.log10(p)
				if logp>=p_thres:
					op_line_sig=line_lst[0]+'\t'+line_lst[1]+'\t'+line_lst[3]+'\t'+line_lst[4]+'\t'+gene_str
					for idx in l_idx_lst:
						gt=line_lst[idx].split(':')[0]
						dp=line_lst[idx].split(':')[1]
						op_line_sig+='\t'+gt+';'+dp

					for idx in s_idx_lst:
						gt=line_lst[idx].split(':')[0]
						dp=line_lst[idx].split(':')[1]	
						op_line_sig+='\t'+gt+';'+dp

					op_line_sig+='\n'
					op2.write(op_line_sig)




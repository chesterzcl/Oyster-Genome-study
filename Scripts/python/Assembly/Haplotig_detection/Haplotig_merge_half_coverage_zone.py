import pandas as pd
import os

os.chdir("/Users/lizicheng/Desktop/Data/dog_genome/oyster/residual_haplotig")

def merge_half_coverage_bed(
    input_bed,
    output_bed,
    merge_distance=4000
    
    min_size=500000
):
    # Load BED
    bed = pd.read_csv(input_bed, sep='\t', header=None, names=['chrom', 'start', 'end'])
    
    merged = []
    for chrom, group in bed.groupby('chrom'):
        group = group.sort_values('start')
        current_start = None
        current_end = None

        for _, row in group.iterrows():
            if current_start is None:
                current_start = row['start']
                current_end = row['end']
            else:
                if row['start'] <= current_end + merge_distance:
                    current_end = max(current_end, row['end'])
                else:
                    merged.append([chrom, current_start, current_end])
                    current_start = row['start']
                    current_end = row['end']
        if current_start is not None:
            merged.append([chrom, current_start, current_end])
    
    # Convert to DataFrame
    merged_df = pd.DataFrame(merged, columns=['chrom', 'start', 'end'])
    merged_df['length'] = merged_df['end'] - merged_df['start']
    
    # Filter for large merged regions
    big_merged = merged_df[merged_df['length'] >= min_size]
    
    # Write output BED
    big_merged[['chrom', 'start', 'end']].to_csv(output_bed, sep='\t', header=False, index=False)
    
    print(f"Merged and filtered BED written to: {output_bed}")



merge_half_coverage_bed("half_coverage_regions.bed", "major_haplotigs.bed", merge_distance=300000, min_size=1000000)


def pick_best_haplotig_per_half_coverage(
    half_coverage_bed,
    self_alignment_tsv,
    output_bed
):
    # Load half-coverage regions
    half_df = pd.read_csv(half_coverage_bed, sep='\t', header=None, names=['chrom', 'start', 'end'])

    # Load filtered self-alignment hits
    align_df = pd.read_csv(self_alignment_tsv, sep='\t')

    # Prepare results
    results = []

    # For each half-coverage zone
    for idx, half in half_df.iterrows():
        chrom = half['chrom']
        h_start = half['start']
        h_end = half['end']

        hits = align_df[
            (align_df['ref_chr'] == chrom) &
            (align_df['qry_chr'] == chrom) &
            ((align_df['ref_end'] > h_start) & (align_df['ref_start'] < h_end)) &
            ((align_df['qry_end'] > h_start) & (align_df['qry_start'] < h_end))
        ].copy()

        hits = hits[hits['alignment_length'] >= 5000].copy()

        if hits.empty:
            continue  # Nothing overlapping this half-coverage region

        # Define alignment span to prioritize biggest
        hits['span'] = hits[['ref_end','ref_start']].max(axis=1) - hits[['ref_end','ref_start']].min(axis=1)
        hits['span_qry'] = hits[['qry_end','qry_start']].max(axis=1) - hits[['qry_end','qry_start']].min(axis=1)
        hits['max_span'] = hits[['span','span_qry']].max(axis=1)

        for _, row in hits.iterrows():
            ref_length = abs(row['ref_end'] - row['ref_start'])
            qry_length = abs(row['qry_end'] - row['qry_start'])
            
            if row['ref_start'] <= row['qry_start']:
                results.append([row['ref_chr'], row['ref_start'], row['ref_end'],row['qry_start'],row['qry_end']])
            else:
                results.append([row['qry_chr'], row['qry_start'], row['qry_end'], row['ref_start'], row['ref_end']])


    # Write final BED
    # Make DataFrame from collected intervals
    result_df = pd.DataFrame(results, columns=['chrom','start','end','start2','end2'])

    # Sort by chrom and start
    result_df = result_df.sort_values(by=['chrom', 'start']).reset_index(drop=True)
    merged_results = []
    merge_distance=3000

    for chrom, group in result_df.groupby('chrom'):
        group = group.sort_values('start')
        current_start = None
        current_end = None

        for _, row in group.iterrows():
            if current_start is None:
                current_start = row['start']
                current_end = row['end']
            elif row['start'] <= current_end + merge_distance:
                current_end = max(current_end, row['end'])
            else:
                merged_results.append([chrom, current_start, current_end])
                current_start = row['start']
                current_end = row['end']
        if current_start is not None:
            merged_results.append([chrom, current_start, current_end])


    result_df.to_csv(output_bed, sep='\t', header=False, index=False)
    print(f"Saved recommended masking BED to: {output_bed}")

    

pick_best_haplotig_per_half_coverage(
    half_coverage_bed="major_haplotigs.bed",
    self_alignment_tsv="coords_large_sv.tsv",
    output_bed="best_haplotigs_to_mask.bed"
)
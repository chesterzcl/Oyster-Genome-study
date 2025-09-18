import pandas as pd
import argparse

def bin_depth(input_file, output_file, bin_size):
    df = pd.read_csv(input_file, sep='\t', header=None, names=['chrom', 'pos', 'depth'])
    df['bin'] = df['pos'] // bin_size
    binned = df.groupby(['chrom', 'bin'])['depth'].mean().reset_index()
    binned['start'] = binned['bin'] * bin_size
    binned['end'] = binned['start'] + bin_size
    binned[['chrom', 'start', 'end', 'depth']].to_csv(output_file, sep='\t', index=False, header=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, help='Input depth file from samtools')
    parser.add_argument('-o', '--output', required=True, help='Output binned file')
    parser.add_argument('-b', '--bin', type=int, default=10000, help='Bin size in base pairs')
    args = parser.parse_args()
    bin_depth(args.input, args.output, args.bin)



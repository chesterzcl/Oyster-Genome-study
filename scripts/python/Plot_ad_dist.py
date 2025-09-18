import matplotlib.pyplot as plt
import numpy as np
import os
import seaborn as sns

# Create an output directory to save the plots
output_dir = "/Users/lizicheng/Desktop/Data/dog_genome/oyster/allele_depth_dist"
os.makedirs(output_dir, exist_ok=True)

# Initialize a dictionary to store allele depths for each sample
# Sample name will be the key, and a list of depths will be the value
sample_allele_depths = {}

cnt=0
col_sample_dict={}
# Read the allele depth data from the file
with open("/Users/lizicheng/Desktop/Data/dog_genome/oyster/20cv_hoh_biSNP_filtered_gt.txt", "r") as file:

    for line in file:
        cnt+=1
        if cnt%100000==0:
            print(cnt)

        # if cnt > 1000000:
        #     break

        line_lst=line.strip().split("\t")

        # Skip the header line if present
        if line_lst[0] == "CHROM":
            for i in range(3,len(line_lst)):
                sample_id=line_lst[i].split('_')[0]
                col_sample_dict[i]=sample_id
            continue

        for i in range(3,len(line_lst)):

            # Split the AD field into depths for all alleles (reference and alternate)
            depths = list(map(int,line_lst[i].split(",")))  # Convert depth strings to integers
            ref_depth=depths[0]
            alt_depth=depths[1]

            # If the sample is not yet in the dictionary, initialize it
            if i not in sample_allele_depths:
                sample_allele_depths[i] = {'reference':{},'alternate':{}}

            if ref_depth not in sample_allele_depths[i]['reference']:
                sample_allele_depths[i]['reference'][ref_depth]=0
            sample_allele_depths[i]['reference'][ref_depth]+=1

            if alt_depth not in sample_allele_depths[i]['alternate']:
                sample_allele_depths[i]['alternate'][alt_depth]=0           
            sample_allele_depths[i]['alternate'][alt_depth]+=1




# Generate and save a histogram for each sample
for sample,depths in sample_allele_depths.items():
    # Get reference and alternate allele depths (dictionaries)

    reference_depths=depths['reference']
    alternate_depths=depths['alternate']
    
    # Sort the dictionaries by allele depth (keys)
    reference_depths = sorted(reference_depths.items())
    alternate_depths = sorted(alternate_depths.items())
    
    # Extract x (allele depths) and y (frequencies) from sorted data
    ref_x, ref_y = zip(*reference_depths)  # Unzip into x and y values for reference
    alt_x, alt_y = zip(*alternate_depths)  # Unzip into x and y values for alternate

    # Filter the x-values (allele depths) to only include those between 5 and 80
    # ref_x_filtered = [x for x in ref_x if 5 <= x <= 100]
    # ref_y_filtered = [y for x, y in zip(ref_x, ref_y) if 5 <= x <= 100]
    
    alt_x_filtered = [x for x in alt_x if 6 <= x <= 100]
    alt_y_filtered = [y for x, y in zip(alt_x, alt_y) if 6 <= x <= 100]
    total=sum(alt_y_filtered)
    for i in range(len(alt_y_filtered)):
        alt_y_filtered[i]=alt_y_filtered[i]/total
    
    # Plot the line plot for the sample
    plt.figure(figsize=(12, 6))
    
    # Plot reference allele depth distribution (filtered)
    # plt.plot(ref_x_filtered, ref_y_filtered, marker='o', linestyle='-', color='blue', label='Reference Allele Depth')

    # Plot alternate allele depth distribution (filtered)
    plt.plot(alt_x_filtered, alt_y_filtered, marker='o', linestyle='-', color='red', label='Variant Allele Depth')

    # Add title and labels
    plt.title(f"Allele Depth Distribution for Sample {col_sample_dict[sample]}", fontsize=16)
    plt.xlabel("Allele Depth", fontsize=14)
    plt.ylabel("Frequency", fontsize=14)

    # Set the x-axis range between 0 and 100
    plt.xlim(6,100)

    # Display grid for easier interpretation
    plt.grid(True)

    # Add a legend to distinguish reference and alternate alleles
    plt.legend()


    # Save the plot as a PNG file
    plt.savefig(os.path.join(output_dir, f"Sample_{col_sample_dict[sample]}_alt_allele_depth_distribution.png"))
    plt.close()

    print(f"Plots saved in {output_dir}/")



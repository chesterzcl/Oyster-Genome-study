import sys
import matplotlib.pyplot as plt



dir="/Users/lizicheng/Desktop/Data/dog_genome/oyster"
data="PB.stat"

cov_list=[]
cnt_list=[]
with open(dir+"/"+data, 'r') as f:
    for line in f:
    # Skip header lines starting with '#'
        if line.startswith('#'):
            continue
        parts = line.strip().split()
        cov=int(parts[0])
        cnt=int(parts[1])
        cov_list.append(cov)
        cnt_list.append(cnt)

# print(cov_dict)

plt.figure(figsize=(10, 6))  # Adjust figure size as needed


plt.plot(cov_list,cnt_list, marker='o', linestyle='-', markersize=1,label="Polished contigs") # Added label


plt.title("HiFi Reads coverage plot")
plt.xlabel("Read depth")
plt.ylabel("Number of bases")
plt.xlim(xmin=0,xmax=300)
plt.grid(True, which="both", linestyle='--', linewidth=0.5)  # Add gridlines
plt.legend()
cutoffs = [5, 42, 96, 97, 175, 411]
labels = ['JUNK', 'HAPLO_MIN', 'HAPLO_MAX', 'DIPLO_MIN', 'DIPLO_MAX', 'HIGHCOV']
colors = ['gray', 'orange', 'orange', 'green', 'green', 'red']

for c, label, color in zip(cutoffs, labels, colors):
    plt.axvline(x=c, color=color, linestyle='--', linewidth=1.5, label=label)
plt.tight_layout() #Adjusts the padding
plt.show()
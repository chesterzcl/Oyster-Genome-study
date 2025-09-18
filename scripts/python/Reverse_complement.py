def reverse_complement(seq):
    complement = str.maketrans('ATCGatcg', 'TAGCtagc')
    return seq.translate(complement)[::-1]

seq1 = "TTAGGG"

print(reverse_complement(seq1))  

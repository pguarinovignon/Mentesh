# A Python program to print all
# combinations of given length
from itertools import combinations
import sys
N= sys.argv[1]

# Get all combinations of [1, 2, 3, etc]
# and length x
comb = combinations(["Iran_GanjDareh_N", "CHG", "PPN", "Barcin_N", "Anatolia_TellKurdu_EC","IRQ_Nemrik9_PPN","TUR_SE_Mardin_PPN"], N)

# Print the obtained combinations
for i in list(comb):
    print(i) 

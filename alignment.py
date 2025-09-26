from Bio import pairwise2
from Bio.pairwise2 import format_alignment

# İki DNA dizisi
seq1 = "ACCGTGTACACGTACGT"
seq2 = "ACGTACGTCAGTACGT"

# Global Alignment : Tüm diziyi hizalar.
print("Global Alignment")
alignments = pairwise2.align.globalxx(seq1, seq2)

for a in alignments:
    print(format_alignment(*a))

#Local Alignment : Sadece benzeyen dizileri hizalar.
print("Local Alignment")
alignments2 = pairwise2.align.localxx(seq1, seq2)

for a in alignments2:
    print(format_alignment(*a))

#Skorlamalı Alignment : Belirli parametrelerle daha biyolojik bir alignment (hizalama) yapar.
alignments3 = pairwise2.align.globalms(seq1, seq2, 2, -1, -0.5, -0.1)
# parametreler: match=2, mismatch=-1, gap open=-0.5, gap extend=-0.1

for a in alignments3:
    print(format_alignment(*a))
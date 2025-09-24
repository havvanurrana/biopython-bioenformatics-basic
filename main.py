#Bu çalışmada temel biyoinformatik çalışması yapılmıştır.

#SeqIO -> Sequence input - output
from Bio import SeqIO

from Bio.SeqUtils import seq3
from Bio.SeqUtils import gc_fraction
from collections import Counter
import matplotlib.pyplot as plt
import pandas as pd
#fasta formatı
for record in SeqIO.parse("brca1.fasta", "fasta"):
    #FASTA formatında "<" ile başlayan ilk kelime genetik kodun ID'sidir.
    print("ID:", record.id)

    print("Uzunluk:", len(record.seq))
    # genetik kodu proteine çevirir
    print(record.seq.translate())
    # tüm genetik kodu yazdırır
    print(record.seq)
    #var olan tüm genetik kodu 5'-3' yönüne göre çevirir
    print(record.seq.reverse_complement())
    #var olan tüm kodu direkt tersine çevirir (TAGC- ATCG şeklinde)
    print(record.seq.complement())
    #var olan tüm kodu RNA'ya çevirir
    print(record.seq.transcribe())
    #proteinleri harflerinin ilk 3'ünü göstermek için
    protein = record.seq.translate(to_stop=True)  # ilk stop kodona kadar çevir
    print(seq3(protein))
    print(len(protein))
    # Aminoasit frekanslarını say
    counts = Counter(protein)

    # Toplam uzunluk
    total = len(protein)

    print("Protein uzunluğu:", total)
    print("Aminoasit kompozisyonu (frekans ve yüzde):\n")

    # DataFrame (tablo) oluştur
    df = pd.DataFrame({
        "Aminoasit": list(counts.keys()),
        "Sayı": list(counts.values())
    })
    df["Yüzde"] = (df["Sayı"] / total) * 100

    for aa, count in counts.items():
        yuzde = (count / total) * 100
        print(f"{aa}: {count} kez, %{yuzde:.2f}")
    # Grafik çiz
    plt.figure(figsize=(10, 5))
    plt.bar(df["Aminoasit"], df["Sayı"])
    plt.xlabel("Aminoasit")
    plt.ylabel("Sayı")
    plt.title("BRCA1 Protein Kompozisyonu")
    plt.show()
    print("GC Oranı:", gc_fraction(record.seq))
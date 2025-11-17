import matplotlib.pyplot as plt
from Bio import SeqIO
import numpy as np

def calc_quality(seq):
    """Simula una 'calidad' como FastQC pero para FASTA."""
    qualities = []
    for base in seq:
        if base in "ACGT":
            qualities.append(40)      # perfecta
        elif base == "N":
            qualities.append(5)       # mala
        else:
            qualities.append(20)      # ambigua
    return qualities

# Leer secuencias
pre = str(next(SeqIO.parse("sequence.fasta", "fasta")).seq)
post = str(next(SeqIO.parse("HPV18_filtrada.fasta", "fasta")).seq)

# Calcular calidades
pre_q = calc_quality(pre)
post_q = calc_quality(post)

# Gráfica
plt.figure(figsize=(14,5))

plt.plot(pre_q, label="Pre-filtrado", alpha=0.8)
plt.plot(post_q, label="Post-filtrado", alpha=0.8)

plt.title("Comparación de calidad simulada por posición (tipo FastQC)")
plt.xlabel("Posición en la secuencia")
plt.ylabel("Calidad (simulada)")
plt.legend()
plt.grid(alpha=0.2)

# Guardar en carpeta QC_output
plt.savefig("QC_output/comparacion_pre_vs_post.png", dpi=300)

print("Gráfico generado en QC_output/comparacion_pre_vs_post.png")

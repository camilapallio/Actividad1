from Bio import SeqIO
import matplotlib.pyplot as plt

# === Cargar secuencia ===
record = SeqIO.read("HPV18_filtrada.fasta", "fasta")
seq = record.seq
seq_str = str(seq)

# === Longitud ===
longitud = len(seq_str)

# === Conteo de nucleótidos ===
conteo = {
    "A": seq_str.count("A"),
    "T": seq_str.count("T"),
    "C": seq_str.count("C"),
    "G": seq_str.count("G")
}

# === %GC ===
gc = (conteo["G"] + conteo["C"]) / longitud * 100

# === Gráfico de composición ===
plt.figure()
plt.bar(conteo.keys(), conteo.values())
plt.xlabel("Nucleótido")
plt.ylabel("Frecuencia")
plt.title("Composición nucleotídica HPV18")
plt.savefig("QC_output/composicion_bases.png")
plt.close()

# === Gráfico %GC ===
plt.figure()
plt.bar(["GC %"], [gc])
plt.ylabel("Porcentaje")
plt.title("Contenido GC HPV18")
plt.savefig("QC_output/gc_content.png")
plt.close()

# === Guardar resumen QC ===
with open("QC_output/QC_resumen.txt", "w") as f:
    f.write(f"Longitud de la secuencia: {longitud}\n")
    f.write(f"Conteo de bases: {conteo}\n")
    f.write(f"GC%: {gc:.2f}%\n")

print("QC generado en carpeta QC_outpu")


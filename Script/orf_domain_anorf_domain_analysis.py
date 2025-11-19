from Bio import SeqIO
from Bio.Seq import Seq
import os


# Archivos de entrada y salida
fasta_file = "GENOMA_HPV18/SECUENCIAS/HPV18_filtrada.fasta"
output_dir = "GENOMA_HPV18/Analisis_secuencia_filtrada"
orf_output_file = os.path.join(output_dir, "orfs_output.txt")
protein_fasta_file = os.path.join(output_dir, "protein_sequences.fasta")

orf_min_len = 100  # longitud mínima de ORFs

# Abrimos archivos de salida
with open(orf_output_file, "w") as orf_out, open(protein_fasta_file, "w") as prot_out:
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = record.seq
        orfs = []
        
        # Buscar ORFs en los 3 marcos de lectura
        for frame in range(3):
            seq_frame = seq[frame:]
            start_codon = "ATG"
            stop_codons = ["TAA", "TAG", "TGA"]
            
            i = 0
            while i < len(seq_frame) - 2:
                codon = seq_frame[i:i+3]
                if codon == start_codon:
                    for j in range(i, len(seq_frame)-2, 3):
                        stop_codon = seq_frame[j:j+3]
                        if stop_codon in stop_codons:
                            orf_seq = seq_frame[i:j+3]
                            if len(orf_seq) >= orf_min_len:
                                start_pos = frame + i + 1
                                end_pos = frame + j + 3
                                orfs.append((start_pos, end_pos, orf_seq))
                            i = j + 3
                            break
                i += 3
        
        # Guardar ORFs en archivo de texto y traducir a proteínas
        orf_out.write(f">{record.id}\n")
        for idx, (start, end, orf_seq) in enumerate(orfs, 1):
            orf_out.write(f"ORF{idx}: {start}-{end}, longitud={len(orf_seq)}\n")
            protein_seq = orf_seq.translate(to_stop=True)
            prot_out.write(f">{record.id}_ORF{idx}\n{protein_seq}\n")

print("✅ Detección de ORFs y traducción completada.")
print(f"Archivos generados en: {output_dir}")
print("Ahora podés usar HMMER o InterProScan sobre protein_sequences.fasta para detectar dominios conservados.")


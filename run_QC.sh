#!/bin/bash
# run_QC.sh
# Uso: ./run_QC.sh <input_fasta> <sequence_id>

set -euo pipefail

IN_FASTA=${1:-sequence.fasta}
SEQ_ID=${2:-AY262282.1}
OUT_DIR="QC_output"

mkdir -p "$OUT_DIR"

echo "Run QC - $(date)" > "$OUT_DIR/QC_stats.txt"
echo "Input fasta: $IN_FASTA" >> "$OUT_DIR/QC_stats.txt"
echo "Sequence ID: $SEQ_ID" >> "$OUT_DIR/QC_stats.txt"
echo "" >> "$OUT_DIR/QC_stats.txt"

samtools faidx "$IN_FASTA"

samtools faidx "$IN_FASTA" "$SEQ_ID" > "$OUT_DIR/HPV18_filtrada.fasta"

seqkit stats "$OUT_DIR/HPV18_filtrada.fasta" >> "$OUT_DIR/QC_stats.txt"

seqkit seq -w 0 "$OUT_DIR/HPV18_filtrada.fasta" > "$OUT_DIR/tmp_one_line.fasta"
grep -v "^>" "$OUT_DIR/tmp_one_line.fasta" | tr -d '\n' > "$OUT_DIR/tmp_seq.txt"

LEN=$(wc -c < "$OUT_DIR/tmp_seq.txt" | tr -d ' ')
NCOUNT=$(grep -o -i "N" "$OUT_DIR/tmp_seq.txt" | wc -l || true)
GCOUNT=$(grep -o -i "G" "$OUT_DIR/tmp_seq.txt" | wc -l || true)
CCOUNT=$(grep -o -i "C" "$OUT_DIR/tmp_seq.txt" | wc -l || true)
GC=$((GCOUNT + CCOUNT))

GC_PERC=$(awk -v g="$GC" -v L="$LEN" 'BEGIN { if (L>0) printf "%.2f", (g/L)*100; else print "NA" }')

echo "" >> "$OUT_DIR/QC_stats.txt"
echo "Length (nt): $LEN" >> "$OUT_DIR/QC_stats.txt"
echo "N count: $NCOUNT" >> "$OUT_DIR/QC_stats.txt"
echo "GC count: $GC (GC% = $GC_PERC)" >> "$OUT_DIR/QC_stats.txt"

cat > "$OUT_DIR/QC_report.md" <<EOF
# QC Report - HPV18 filtered sequence
Fecha: $(date)

## Archivo analizado
- Input: $IN_FASTA
- Secuencia filtrada: HPV18_filtrada.fasta (ID: $SEQ_ID)

## Resultados QC Rápidos
- Longitud: $LEN nt
- Conteo de N: $NCOUNT
- %GC: $GC_PERC%

## Archivos generados
- HPV18_filtrada.fasta
- QC_stats.txt
- QC_report.md

EOF

echo "QC completado. Salida en $OUT_DIR/"

# ============================================================
# Extracción de promotores + análisis de motivos con MEME
# ============================================================

# --- 1. Librerías necesarias ---
# BiocManager::install(c("rtracklayer", "Biostrings", "memes", "universalmotif"))

library(rtracklayer)
library(Biostrings)
library(memes)
library(universalmotif)

# --- 2. Rutas de archivos (ajusta si es necesario) ---
genome_fasta   <- "GCF_036785885.1_Coffea_Arabica_ET-39_HiFi_genomic.fna"
annotation_gff <- "genomic.gff"

file.exists(genome_fasta)     # debe ser TRUE
file.exists(annotation_gff)   # debe ser TRUE

# --- 3. Genes de interés ---
genes_interes <- c(
  "LOC113695805", "LOC140009270", "LOC113690379", "LOC113690378",
  "LOC113699104", "LOC113691088", "LOC113691272", "LOC113688572"
)

# --- 4. Importar el GFF y quedarte solo con filas tipo "gene" ---
gff <- import(annotation_gff, format = "gff3")
genes_gff <- gff[gff$type == "gene"]

# Revisa qué columna trae el ID del gen (corre esto UNA VEZ para confirmar)
head(mcols(genes_gff))

# --- 5. Filtrar tus genes de interés ---
# Ajusta esta línea según lo que veas en el head() de arriba:

# Opción A: columna "gene" tiene el ID limpio (ej. "LOC113695805")
genes_gr_filtrados <- genes_gff[genes_gff$gene %in% genes_interes]

# Opción B (usar solo si la A no encuentra nada): ID viene como "gene-LOC113695805"
# genes_gr_filtrados <- genes_gff[sub("gene-", "", genes_gff$ID) %in% genes_interes]

cat("Genes encontrados:", length(genes_gr_filtrados), "de", length(genes_interes), "\n")
print(genes_gr_filtrados$gene)   # cambia a $ID si usaste la Opción B

# --- 6. Definir la región promotora (respeta hebra +/- automáticamente) ---
upstream_bp   <- 1000
downstream_bp <- 200

promotores_gr <- promoters(genes_gr_filtrados,
                           upstream = upstream_bp,
                           downstream = downstream_bp)
library(Rsamtools)

# --- 7. Cargar el genoma como FaFile indexado ---
# Esto crea un índice .fai la primera vez (puede tardar un poco con genomas grandes)
if (!file.exists(paste0(genome_fasta, ".fai"))) {
  indexFa(genome_fasta)
}

fa <- FaFile(genome_fasta)

# --- Asegurarse de que los seqnames coincidan entre el GFF y el genoma ---
seqlevels(promotores_gr)          # nombres de cromosomas/scaffolds en tu GRanges
seqnames(scanFaIndex(fa))         # nombres de cromosomas/scaffolds en el FASTA indexado

# Si no coinciden exactamente, hay que ajustar seqlevels(promotores_gr) antes de continuar

# --- Ajustar seqinfo con las longitudes reales del genoma indexado ---
fa_seqinfo <- seqinfo(fa)
seqlevels(promotores_gr) <- seqlevels(promotores_gr)[seqlevels(promotores_gr) %in% seqnames(fa_seqinfo)]
seqlengths(promotores_gr) <- seqlengths(fa_seqinfo)[seqlevels(promotores_gr)]
promotores_gr <- trim(promotores_gr)

# --- Extraer las secuencias promotoras ---
promotores_seq <- getSeq(fa, promotores_gr)
names(promotores_seq) <- genes_gr_filtrados$gene   # o $ID, según cuál hayas usado

promotores_seq
writeXStringSet(promotores_seq, filepath = "promotores_8genes.fasta")

options(meme_bin = "/opt/miniconda3/envs/meme_env/bin/")
check_meme_install()
# --- 8. Correr MEME de novo ---
check_meme_install()

resultados_meme <- runMeme(
  input   = promotores_seq,
  db      = NULL,
  nmotifs = 3,
  minw    = 6,
  maxw    = 15,
  mod     = "zoops",
  outdir  = "meme_output_8genes",
  parse_genomic_coord = FALSE
)

print(resultados_meme)

library(universalmotif)

# Convertir de data.frame a lista de objetos universalmotif
motivos_um <- to_list(resultados_meme)

# Ahora sí visualizar
view_motifs(motivos_um)

jaspar_db <- "JASPAR_plants.meme"
file.exists(jaspar_db)   # debe ser TRUE

resultados_tomtom <- runTomTom(
  input    = resultados_meme,
  database = jaspar_db
)

print(resultados_tomtom)

library(dplyr)
install.packages("dplyr")
library(dplyr)
# Ver los top matches con su score de similitud
resultados_tomtom %>%
  select(name, best_match_name, best_match_altname, best_match_pval, best_match_qval)

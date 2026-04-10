import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import matplotlib.pyplot as plt
from collections import defaultdict

def load_snps(path):
    snps = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if (
                line.startswith("/") or
                line.startswith("NUCMER") or
                line.startswith("[") or
                line.startswith("=") or
                line == ""
            ):
                continue

            parts = line.split()
            if len(parts) < 3:
                continue

            try:
                P1 = int(parts[0])
            except ValueError:
                continue  # если вдруг встретится мусор

            ref_nt = parts[1]
            alt_nt = parts[2]

            snps.append((P1, ref_nt, alt_nt))

    df = pd.DataFrame(snps, columns=["pos", "ref", "alt"])
    return df

snp_df = load_snps("ecoli_snps_for_report.txt")

gb_record = SeqIO.read("ecoli_k12.gb", "genbank")
genome_seq = gb_record.seq

genes = []
for feat in gb_record.features:
    if feat.type == "CDS":
        start = int(feat.location.start)
        end = int(feat.location.end)
        strand = feat.location.strand
        gene = feat.qualifiers.get("gene", ["unknown"])[0]
        product = feat.qualifiers.get("product", ["unknown"])[0]
        genes.append((start, end, strand, gene, product))

genes_df = pd.DataFrame(genes, columns=["start", "end", "strand", "gene", "product"])

def find_gene(pos):
    hits = genes_df[(genes_df.start <= pos) & (genes_df.end >= pos)]
    if len(hits) == 0:
        return None
    return hits.iloc[0]

def classify_snp(gene_row, pos, ref, alt):
    if ref == "." or alt == ".":
        return "indel", None

    if gene_row is None:
        return "intergenic", None

    start, end, strand, gene, product = gene_row

    cds_pos = pos - start
    if strand == -1:
        cds_pos = end - pos

    codon_index = cds_pos // 3
    codon_pos = cds_pos % 3

    cds_seq = genome_seq[start:end]
    if strand == -1:
        cds_seq = cds_seq.reverse_complement()

    codon = list(str(cds_seq[codon_index*3:(codon_index+1)*3]))

    if len(codon) != 3:
        return "partial", None

   codon[codon_pos] = alt
    new_codon = "".join(codon)

    try:
        aa_old = str(Seq("".join(codon)).translate())
        aa_new = str(Seq(new_codon).translate())
    except:
        return "invalid", None

    if aa_old == aa_new:
        snp_type = "synonymous"
    else:
        snp_type = "nonsynonymous"

    aa_change = f"p.{aa_old}{codon_index+1}{aa_new}"

    return snp_type, aa_change

annot_gene = []
annot_type = []
annot_aa = []
annot_product = []

for i, row in snp_df.iterrows():
    pos, ref, alt = row["pos"], row["ref"], row["alt"]

    gene_row = find_gene(pos)
    if gene_row is None:
        annot_gene.append("intergenic")
        annot_type.append("intergenic")
        annot_aa.append(None)
        annot_product.append(None)
        continue

    gene = gene_row.gene
    product = gene_row.product

    snp_type, aa_change = classify_snp(gene_row, pos, ref, alt)

    annot_gene.append(gene)
    annot_type.append(snp_type)
    annot_aa.append(aa_change)
    annot_product.append(product)

snp_df["gene"] = annot_gene
snp_df["type"] = annot_type
snp_df["aa_change"] = annot_aa
snp_df["product"] = annot_product

def functional_effect(row):
    if row["type"] == "intergenic":
        return "likely regulatory"
    if row["type"] == "synonymous":
        return "likely neutral"
    if row["type"] == "nonsynonymous":
        return "possible functional impact"
    return "unknown"

snp_df["effect"] = snp_df.apply(functional_effect, axis=1)

snp_df.to_csv("annotated_snps.csv", index=False)

plt.figure(figsize=(14, 4))
plt.hist(snp_df["pos"], bins=200, color="steelblue")
plt.title("SNP density along genome")
plt.xlabel("Genome position")
plt.ylabel("Count")
plt.tight_layout()
plt.savefig("snp_density_final.png", dpi=300)

plt.figure(figsize=(6, 4))
snp_df["type"].value_counts().plot(kind="bar", color="darkred")
plt.title("Distribution of SNP types")
plt.ylabel("Count")
plt.tight_layout()
plt.savefig("snp_types.png", dpi=300)

top10 = (
    snp_df[snp_df.gene != "intergenic"]
    .groupby("gene")
    .size()
    .sort_values(ascending=False)
    .head(10)
)

top10.to_csv("top10_genes.csv")
print("Top-10 genes saved to top10_genes.csv")

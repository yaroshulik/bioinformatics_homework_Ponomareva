from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd

# ========== ЗАГРУЗКА ДАННЫХ ==========
print("Загрузка K-12...")
k12 = SeqIO.read("ecoli_k12.gb", "genbank")
print("Загрузка O157:H7...")
o157 = SeqIO.read("ecoli_o157.gb", "genbank")

def extract_proteins(gb_record):
    """Извлекает все белки из GenBank записи"""
    proteins = {}
    for feat in gb_record.features:
        if feat.type == "CDS":
            # Получаем ген и продукт
            gene = feat.qualifiers.get("gene", ["unknown"])[0]
            product = feat.qualifiers.get("product", ["unknown"])[0]
            
            # Получаем белковую последовательность
            if "translation" in feat.qualifiers:
                protein_seq = feat.qualifiers["translation"][0]
                proteins[gene] = {
                    "product": product,
                    "sequence": protein_seq,
                    "location": f"{int(feat.location.start)}-{int(feat.location.end)}"
                }
    return proteins

# Извлекаем белки из обоих геномов
k12_proteins = extract_proteins(k12)
o157_proteins = extract_proteins(o157)

print(f"Найдено белков в K-12: {len(k12_proteins)}")
print(f"Найдено белков в O157:H7: {len(o157_proteins)}")

# ========== АНАЛИЗ РАЗЛИЧИЙ ==========
results = []

# 1. Уникальные гены O157 (нет в K-12)
o157_unique = set(o157_proteins.keys()) - set(k12_proteins.keys())
print(f"\nУникальных генов в O157:H7: {len(o157_unique)}")

for gene in o157_unique:
    results.append({
        "gene": gene,
        "type": "UNIQUE_TO_O157",
        "product": o157_proteins[gene]["product"],
        "location_o157": o157_proteins[gene]["location"],
        "aa_changes": "GENE_ABSENT_IN_K12"
    })

# 2. Уникальные гены K-12 (нет в O157)
k12_unique = set(k12_proteins.keys()) - set(o157_proteins.keys())
print(f"Уникальных генов в K-12: {len(k12_unique)}")

for gene in k12_unique:
    results.append({
        "gene": gene,
        "type": "UNIQUE_TO_K12",
        "product": k12_proteins[gene]["product"],
        "location_o157": "N/A",
        "aa_changes": "GENE_ABSENT_IN_O157"
    })

# 3. Общие гены с различиями
common_genes = set(k12_proteins.keys()) & set(o157_proteins.keys())
print(f"Общих генов: {len(common_genes)}")

nonsyn_count = 0
for gene in common_genes:
    k12_seq = k12_proteins[gene]["sequence"]
    o157_seq = o157_proteins[gene]["sequence"]
    
    if k12_seq != o157_seq:
        # Находим позиции различий
        changes = []
        for i, (aa_k12, aa_o157) in enumerate(zip(k12_seq, o157_seq)):
            if aa_k12 != aa_o157:
                changes.append(f"p.{aa_k12}{i+1}{aa_o157}")
        
        # Определяем тип изменений
        if len(k12_seq) != len(o157_seq):
            snp_type = "LENGTH_DIFFERENCE"
        else:
            snp_type = "NONSYNONYMOUS"
            nonsyn_count += 1
        
        results.append({
            "gene": gene,
            "type": snp_type,
            "product": o157_proteins[gene]["product"],
            "location_o157": o157_proteins[gene]["location"],
            "aa_changes": ", ".join(changes) if changes else "LENGTH_MISMATCH"
        })

print(f"Общих генов с несинонимичными заменами: {nonsyn_count}")

# ========== СОХРАНЕНИЕ РЕЗУЛЬТАТОВ ==========
df = pd.DataFrame(results)

# Сохраняем полный отчет
df.to_csv("k12_vs_o157_comparison.csv", index=False)

# Отдельно сохраняем уникальные гены O157 (важные для патогенности!)
o157_unique_df = df[df["type"] == "UNIQUE_TO_O157"]
o157_unique_df.to_csv("o157_unique_genes.csv", index=False)

# Отдельно сохраняем гены с мутациями
mutated_df = df[df["type"] == "NONSYNONYMOUS"]
mutated_df.to_csv("o157_mutations.csv", index=False)

# ========== ВЫВОД СТАТИСТИКИ ==========
print("\n" + "="*50)
print("РЕЗУЛЬТАТЫ СРАВНЕНИЯ K-12 vs O157:H7")
print("="*50)
print(f"Уникальные гены O157:H7 (патогенность, токсины): {len(o157_unique)}")
print(f"Гены с несинонимичными заменами: {nonsyn_count}")
print(f"Гены, отсутствующие в O157:H7: {len(k12_unique)}")
print("\nТоп-10 уникальных генов O157:H7:")
for gene in list(o157_unique)[:10]:
    print(f"  - {gene}: {o157_proteins[gene]['product'][:60]}...")

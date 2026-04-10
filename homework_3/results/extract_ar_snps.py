import pandas as pd

#ЗАГРУЗКА ДАННЫХ
print("Загрузка файлов...")

#загружаем аннотированные SNP
df = pd.read_csv("annotated_snps.csv")
print(f"  Всего SNP: {len(df)}")

#загружаем список генов устойчивости (пропускаем пустые строки и комментарии)
ar_genes = []
with open("ar_genes.txt") as f:
    for line in f:
        line = line.strip()
        if line and not line.startswith("#"):
            ar_genes.append(line.lower())  #приводим к нижнему регистру для поиска

print(f"  Генов в списке: {len(ar_genes)}")

#ПОИСК SNP В ГЕНАХ УСТОЙЧИВОСТИ
print("\nПоиск SNP в генах устойчивости...")

#приводим колонку с генами к нижнему регистру для сравнения
df["gene_lower"] = df["gene"].str.lower()

#фильтруем: ген должен быть в списке и не NaN
ar_snps = df[df["gene_lower"].isin(ar_genes) & df["gene"].notna()]

#удаляем временную колонку
ar_snps = ar_snps.drop("gene_lower", axis=1)

print(f"  Найдено SNP: {len(ar_snps)}")

#РАЗДЕЛЕНИЕ ПО ТИПАМ
ar_nonsyn = ar_snps[ar_snps["type"] == "nonsynonymous"]
ar_syn = ar_snps[ar_snps["type"] == "synonymous"]
ar_indel = ar_snps[ar_snps["type"] == "indel"]

print(f"    - несинонимичных: {len(ar_nonsyn)}")
print(f"    - синонимичных: {len(ar_syn)}")
print(f"    - инделов: {len(ar_indel)}")

#СОХРАНЕНИЕ РЕЗУЛЬТАТОВ
print("\nСохранение результатов...")

ar_snps.to_csv("ar_all_snps.csv", index=False)
ar_nonsyn.to_csv("ar_nonsyn_snps.csv", index=False)

print("  ar_all_snps.csv - все SNP в генах устойчивости")
print("  ar_nonsyn_snps.csv - только несинонимичные замены")

#ВЫВОД В КОНСОЛЬ
print("\n" + "="*50)
if len(ar_nonsyn) > 0:
    print("НЕСИНОНИМИЧНЫЕ ЗАМЕНЫ В ГЕНАХ УСТОЙЧИВОСТИ:")
    print("="*50)
    
    cols = ["pos", "gene", "aa_change", "product"]
    for _, row in ar_nonsyn.iterrows():
        print(f"\nГен: {row['gene']}")
        print(f"  Позиция в геноме: {row['pos']}")
        print(f"  Замена: {row['aa_change']}")
        print(f"  Продукт: {row['product']}")
else:
    print("РЕЗУЛЬТАТ:")
    print("  Несинонимичных замен в генах устойчивости не обнаружено")
    print("  Штамм вероятно сохраняет чувствительность к антибиотикам")
print("="*50)

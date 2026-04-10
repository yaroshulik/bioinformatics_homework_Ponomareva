import subprocess
import matplotlib.pyplot as plt
import numpy as np

def get_alignment_lengths(delta_file):
    """Извлекает длины выравниваний из .delta через show-coords"""
    cmd = f"show-coords -r {delta_file}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    lines = result.stdout.strip().split('\n')
    
    lengths = []
    for line in lines[5:]:
        if line.strip():
            parts = line.split()
            if len(parts) >= 8:
                lengths.append(int(parts[7]))
    return lengths

print("Загрузка данных...")
lengths_f1 = get_alignment_lengths("ecoli_full_filt1.delta")
lengths_f2 = get_alignment_lengths("ecoli_full_filt2.delta")

print(f"Filter 1: найдено {len(lengths_f1)} выравниваний")
print(f"Filter 2: найдено {len(lengths_f2)} выравниваний")

fig, ax = plt.subplots(figsize=(12, 6))

ax.hist(lengths_f1, bins=50, alpha=0.6, label=f"Filter 1 (id>90%, 1-to-1)\nn={len(lengths_f1)}", 
        color='#66c2a5', edgecolor='black', linewidth=0.5)

ax.hist(lengths_f2, bins=50, alpha=0.6, label=f"Filter 2 (len>5kb, 1-to-1)\nn={len(lengths_f2)}", 
        color='#fc8d62', edgecolor='black', linewidth=0.5)

ax.set_xlabel("Alignment length (bp)", fontsize=12)
ax.set_ylabel("Frequency", fontsize=12)
ax.set_title("Distribution of alignment lengths: E. coli K-12 vs O157:H7", fontsize=14, fontweight='bold')
ax.set_xscale('log')
ax.legend(loc='upper right', fontsize=10)
ax.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig("alignment_lengths_histograms_1.png", dpi=300)
print(f"\n Сохранён график: alignment_lengths_histograms.png")

# Сохраняем статистику в CSV
import csv
with open('alignment_stats.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['Filter', 'Count', 'Min_bp', 'Q1_bp', 'Median_bp', 'Q3_bp', 'Max_bp'])
    writer.writerow(['Filter_1_id90', len(lengths_f1), np.min(lengths_f1), np.percentile(lengths_f1, 25), 
                     np.median(lengths_f1), np.percentile(lengths_f1, 75), np.max(lengths_f1)])
    writer.writerow(['Filter_2_len5kb', len(lengths_f2), np.min(lengths_f2), np.percentile(lengths_f2, 25), 
                     np.median(lengths_f2), np.percentile(lengths_f2, 75), np.max(lengths_f2)])
print(" Сохранена статистика: alignment_stats.csv")

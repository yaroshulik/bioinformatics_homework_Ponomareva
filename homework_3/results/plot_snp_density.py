import subprocess
import matplotlib.pyplot as plt
import numpy as np
import os

def extract_snp_positions(delta_file, filter_name):
    """Извлекает позиции SNP из delta-файла"""
    print(f"Извлечение SNP из {filter_name}...")
    cmd = f"show-snps -r {delta_file} > snps_{filter_name}.txt"
    subprocess.run(cmd, shell=True)
    
    snp_positions = []
    with open(f"snps_{filter_name}.txt", "r") as f:
        for line in f:
            if line.strip() and not line.startswith("NUCMER") and not line.startswith("="):
                parts = line.split()
                if len(parts) >= 1 and parts[0].isdigit():
                    snp_positions.append(int(parts[0]))
    
    print(f"    Найдено SNP: {len(snp_positions)}")
    return snp_positions

def calculate_snp_density(snp_positions, window_kb=1):
    """Рассчитывает плотность SNP на kb (окнами без сглаживания)"""
    if len(snp_positions) == 0:
        return [], [], 0
    
    window_bp = window_kb * 1000
    genome_len = max(snp_positions)
    num_windows = (genome_len // window_bp) + 1
    
    density = []
    window_centers_mbp = []
    
    for i in range(num_windows):
        start = i * window_bp
        end = (i+1) * window_bp
        count = sum(1 for p in snp_positions if start <= p < end)
        density.append(count)
        window_centers_mbp.append((start + end) / 2 / 1e6)
    
    mean_density = np.mean(density)
    return window_centers_mbp, density, mean_density

fig, axes = plt.subplots(3, 1, figsize=(14, 10))
filters = [
    {"file": "ecoli_full_filt1.delta", "name": "Filter 1 (id>90%, 1-to-1)", "color": "steelblue", "mean_color": "darkblue"},
    {"file": "ecoli_full_filt2.delta", "name": "Filter 2 (len>5kb, 1-to-1)", "color": "coral", "mean_color": "darkred"},
    {"file": "ecoli_full_multi.delta", "name": "Filter 3 (multi, id>85%)", "color": "seagreen", "mean_color": "darkgreen"}
]

for idx, f in enumerate(filters):
    print(f"\nОбработка {f['name']}...")
    
    # Извлекаем SNP
    snp_positions = extract_snp_positions(f["file"], f"filter{idx+1}")
    
    if len(snp_positions) == 0:
        print(f"  Предупреждение: нет SNP для {f['name']}")
        continue
    
    # Рассчитываем плотность SNP (окно 10 kb)
    positions_mbp, density_per_kb, mean_density = calculate_snp_density(snp_positions, window_kb=10)
    
    # Рисуем график
    ax = axes[idx]
    
    # Столбчатая диаграмма (гистограмма) плотности
    ax.bar(positions_mbp, density_per_kb, width=0.01, color=f['color'], 
           alpha=0.7, edgecolor='black', linewidth=0.3, label="SNP density")
    
    # Пунктирная линия среднего
    ax.axhline(y=mean_density, color=f['mean_color'], linestyle='--', 
               linewidth=1.5, alpha=0.8, label=f"Mean = {mean_density:.2f} SNPs/kb")
    
    # Настройки графика
    ax.set_xlabel("Position (Mbp)", fontsize=11)
    ax.set_ylabel("SNP density (per kb)", fontsize=11)
    ax.set_title(f"{f['name']}", fontsize=12, fontweight='bold')
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(True, alpha=0.3, axis='y')
    
    # Добавляем информацию
    ax.text(0.02, 0.95, f"Total SNPs: {len(snp_positions)}\nWindow: 10 kb\nMean density: {mean_density:.2f} SNPs/kb", 
            transform=ax.transAxes, fontsize=8, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
    
    # Устанавливаем пределы
    if positions_mbp:
        ax.set_xlim(0, max(positions_mbp))
    if density_per_kb:
        ax.set_ylim(0, max(density_per_kb) * 1.2 if max(density_per_kb) > 0 else 5)

plt.suptitle("SNP Density Along the Genome: E. coli K-12 vs O157:H7\n(Window size = 10 kb)", 
             fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig("snp_density_correct.png", dpi=300, bbox_inches='tight')
print("\n Сохранён график: snp_density_correct.png")

# Также сохраняем статистику
with open("snp_density_stats_correct.csv", "w") as csvfile:
    csvfile.write("Filter,Total_SNPs,Mean_SNP_density_per_kb,Window_kb\n")
    for idx, f in enumerate(filters):
        snp_positions = extract_snp_positions(f["file"], f"filter{idx+1}")
        if len(snp_positions) > 0:
            _, _, mean_density = calculate_snp_density(snp_positions, window_kb=10)
            csvfile.write(f"{f['name']},{len(snp_positions)},{mean_density:.2f},10\n")
print(" Сохранена статистика: snp_density_stats_correct.csv")

# Очистка
for i in range(1, 4):
    if os.path.exists(f"snps_filter{i}.txt"):
        os.remove(f"snps_filter{i}.txt")
print(" Временные файлы удалены")

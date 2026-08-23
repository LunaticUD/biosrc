#!/usr/bin/env python3
# plot_cluster_distribution.py
import matplotlib.pyplot as plt
import numpy as np
# 读取聚类大小分布
sizes = []
with open('cluster_size_distribution.txt', 'r') as f:
	for line in f:
		count, size = map(int, line.strip().split())
		sizes.extend([size] * count)
# 绘制分布图
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
# 直方图
axes[0].hist(sizes, bins=50, alpha=0.7, color='skyblue', edgecolor='black')
axes[0].set_xlabel('Cluster Size')
axes[0].set_ylabel('Frequency')
axes[0].set_title('Cluster Size Distribution')
axes[0].grid(True, alpha=0.3)
# 累积分布
sorted_sizes = np.sort(sizes)
cumulative = np.arange(1, len(sorted_sizes)+1) / len(sorted_sizes)
axes[1].plot(sorted_sizes, cumulative, 'r-', linewidth=2)
axes[1].set_xlabel('Cluster Size')
axes[1].set_ylabel('Cumulative Proportion')
axes[1].set_title('Cumulative Distribution')
axes[1].grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('cluster_distribution.png', dpi=300, bbox_inches='tight')
plt.close()
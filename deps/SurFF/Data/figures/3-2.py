import numpy as np
import matplotlib.pyplot as plt

# Metrics and scores
metrics = ["ID", "OOD"]
dft_scores = np.array([10919, 7576])
surff_scores = np.array([0.23, 0.27])

# Setting the positions and width for the bars
pos = np.arange(len(metrics))
bar_width = 0.25

# fontsize 18
plt.rcParams.update({'font.size': 28,})
plt.figure(figsize=(7, 6))

# Plotting ID scores
plt.bar(pos - bar_width / 2, dft_scores, bar_width, label='DFT', color='darkkhaki', hatch='/')
# add text above bar
for i, v in enumerate(dft_scores):
    plt.text(i - bar_width/2, v + 0.05, str(round(v, 0)), color='black', ha='center',
             fontdict={ 'fontsize': 24})

# Plotting OOD scores
plt.bar(pos + bar_width / 2, surff_scores, bar_width, label='SurFF', color='khaki', hatch='|')
# add text above bar
for i, v in enumerate(surff_scores):
    plt.text(i + bar_width/2, v + 0.05, str(round(v, 2)), color='black', ha='center',
             fontdict={ 'fontsize': 24})

# Adding the aesthetics
# plt.xlabel('Metric')
plt.ylabel('Time (hr)')
plt.title('Computation Time')
plt.xticks(pos, metrics)
plt.yscale('log')
plt.ylim(0.01, 1e5)
plt.legend(fontsize=24, loc='center', framealpha=0.5, bbox_to_anchor=(0.42, 0.7))

# Show the plot
plt.tight_layout()
plt.show()

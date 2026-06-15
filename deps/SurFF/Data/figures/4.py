import numpy as np
import matplotlib.pyplot as plt

# Performance metrics
metrics = ['SP Energy↓\n(eV)',
           'SP Force↓\n(eV/$\AA$)',
           'Relaxed Energy↓\n(eV/$\AA^2$)',
           'Top3 Rate↑\n(%)',
           'Top5 Rate↑\n(%)']
model_a_scores = np.array([2.34, 0.194, 10.5,  66.7, 81.0,])
model_b_scores = np.array([1.29, 0.090, 6.8, 70.0,	82.0])

mean_ab = np.mean([model_a_scores, model_b_scores], axis=0)
mean_ab[-1] = mean_ab[-1]
a_norm = model_a_scores / mean_ab
b_norm = model_b_scores / mean_ab

# Setting the positions and width for the bars
pos = np.arange(len(metrics))
bar_width = 0.15

# fontsize 20
plt.rcParams.update({'font.size': 18,})
plt.figure(figsize=(9, 8))

# Plotting Model A scores
plt.barh(pos - bar_width/1.4, a_norm, bar_width, label='Base', color='cadetblue')
# add text above bar
for i, v in enumerate(a_norm):
    plt.text(v + 0.02, i - bar_width/1.4, str(round(model_a_scores[i], 2)), color='cadetblue', va='center',
             fontdict={'fontweight': 'bold'})

# Plotting Model B scores
plt.barh(pos + bar_width/1.4, b_norm, bar_width, label='Finetune', color = 'mediumseagreen')
# add text above bar
for i, v in enumerate(b_norm):
    plt.text(v + 0.02, i + bar_width/1.4, str(round(model_b_scores[i], 2)), color='mediumseagreen', va='center',
             fontdict={'fontweight': 'bold'})

# Adding the aesthetics
plt.xlabel('Metric Values')
# remove x ticks
plt.xticks([])
plt.ylabel('Metric')
plt.title('Model Performance Comparison')
plt.yticks(pos, metrics)
plt.xlim(0.3, 1.6)
plt.legend(fontsize=16, loc='upper right', framealpha=0.5)

# Show the plot
plt.tight_layout()
plt.show()

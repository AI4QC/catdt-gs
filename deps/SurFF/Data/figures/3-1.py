import numpy as np
import matplotlib.pyplot as plt

# Metrics and scores
metrics = ["Overall", "Low", "Medium", "High", "Top-3", "Top-5"]
id_scores = np.array([0.719, 0.695, 0.692, 0.771, 0.758, 0.8001])
ood_scores = np.array([0.669, 0.666, 0.605, 0.744, 0.666, 0.8101])

# Setting the positions and width for the bars
pos = np.arange(len(metrics))
bar_width = 0.35

# fontsize 18
plt.rcParams.update({'font.size': 28,})
plt.figure(figsize=(14, 6))

# Plotting ID scores
plt.bar(pos - bar_width/2, id_scores, bar_width, label='ID', color='cadetblue',hatch='/')
# add text above bar
for i, v in enumerate(id_scores):
    plt.text(i - bar_width/2 -0.08, v+0.005, f"{v:.3f}", color='black', ha='center',
             fontdict={'fontsize': 24})

# Plotting OOD scores
plt.bar(pos + bar_width/2, ood_scores, bar_width, label='OOD', color='mediumseagreen',hatch='|')
# add text above bar
for i, v in enumerate(ood_scores):
    plt.text(i + bar_width/2 +0.08, v + 0.005, f"{v:.3f}", color='black', ha='center',
             fontdict={'fontsize': 24})

# Adding the aesthetics
# plt.xlabel('Metric')
plt.ylabel('Accuracy')
plt.title('Surface Exposure Accuracy')
plt.xticks(pos, metrics)
plt.ylim(0.5, 0.85)
# place legend at (0.1, 0.4)
plt.legend(fontsize=24,bbox_to_anchor=(0, 1), loc='upper left', framealpha=0.5)

# Show the plot
plt.tight_layout()
plt.show()

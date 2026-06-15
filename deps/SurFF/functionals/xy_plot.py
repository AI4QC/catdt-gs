# plot error in df
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import gaussian_kde


def xy_plot(df, save_pth=None):
    # Example data
    true_values = df["surface_energy_true"].values * 1000
    predicted_values = df["surface_energy_pred"].values * 1000
    errors = predicted_values - true_values
    mae = np.mean(np.abs(errors))

    # Calculate the point density for the scatter plot
    xy = np.vstack([true_values, predicted_values])
    z = gaussian_kde(xy)(xy)


    # set font size
    plt.rcParams.update({'font.size': 14})

    # Create the main scatter plot, marginal plots, and residual plot
    fig = plt.figure(figsize=(6, 6))
    scatter_width = 0.65
    scatter_height = 0.65
    margin_size = 0.05
    residual_height = 0.1
    ax_scatter = plt.axes([0.2, 0.3, scatter_width, scatter_height])
    ax_marg_x = plt.axes([0.2, 0.95, scatter_width, margin_size])
    ax_marg_y = plt.axes([0.85, 0.3, margin_size, scatter_height])
    ax_residual = plt.axes([0.2, 0.28 - residual_height, scatter_width, residual_height], sharex=ax_scatter)

    # Main scatter plot
    ax_scatter.scatter(true_values, predicted_values, c=z, cmap="viridis", s=15, alpha=0.5)
    ax_scatter.set_ylabel('Predicted Values\n (meV/A$^2$)', fontsize=18, fontweight='bold')
    ax_scatter.tick_params(axis='y', labelsize=14)

    # text mae
    ax_scatter.text(0.95, 0.05, f"MAE: {mae:.1f} meV/A$^2$", transform=ax_scatter.transAxes, fontsize=18, fontweight='bold',
                    verticalalignment='bottom', horizontalalignment='right')

    # Add diagonal line
    ax_scatter.plot([min(true_values), max(predicted_values)], [min(true_values), max(predicted_values)], '--'
                    , c='gray', alpha=0.7)

    # KDE plot for x margin
    sns.kdeplot(true_values, ax=ax_marg_x, fill=True)
    ax_marg_x.set(xticks=[], yticks=[], xlim=ax_scatter.get_xlim(), xlabel='', ylabel='')

    # KDE plot for y margin
    sns.kdeplot(y=predicted_values, ax=ax_marg_y, fill=True)
    ax_marg_y.set(xticks=[], yticks=[], ylim=ax_scatter.get_ylim(), xlabel='', ylabel='')

    # Residual plot
    ax_residual.scatter(true_values, errors, alpha=0.2, s=15)
    ax_residual.axhline(0, color='gray', linestyle='--', alpha=0.7)
    ax_residual.axhline(0.01, color='gray', linestyle='--', alpha=0.7)
    ax_residual.axhline(-0.01, color='gray', linestyle='--', alpha=0.7)
    ax_residual.set_xlabel('True Values (meV/A$^2$)', fontsize=18, fontweight='bold')
    ax_residual.set_ylabel('Errors\n (meV/A$^2$)', fontsize=18, fontweight='bold')
    ax_residual.set_yticks([-10, 0, 10])
    ax_residual.set_ylim(-20, 20)
    ax_residual.tick_params(axis='y', labelsize=14)
    ax_residual.tick_params(axis='x', labelsize=14)

    # Hide unnecessary spines and ticks
    for ax in [ax_marg_x, ax_marg_y]:
        for spine in ax.spines.values():
            spine.set_visible(False)

    ax_scatter.spines['top'].set_visible(False)
    ax_scatter.spines['right'].set_visible(False)

    # Synchronize the x-axis limits and hide x-ticks of the scatter plot
    plt.setp(ax_scatter.get_xticklabels(), visible=False)
    plt.tight_layout()

    plt.show()

    if save_pth is not None:
        fig.savefig(save_pth, dpi=600)

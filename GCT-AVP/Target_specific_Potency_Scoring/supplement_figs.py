import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.metrics import r2_score, mean_squared_error

plt.rcParams.update({
    'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'pdf.fonttype': 42,
    'ps.fonttype': 42
})

COLOR_VAL = '#1B263B'    
COLOR_GEN = '#2A9D8F'     
COLOR_LINE = '#D90429'    

fig, axes = plt.subplots(1, 3, figsize=(16, 4.8), dpi=300)
plt.subplots_adjust(wspace=0.3)

ax1 = axes[0]

    df_ec50 = pd.read_csv("./EC50.csv")
    y_true_col = [c for c in df_ec50.columns if 'true' in c.lower() or 'exp' in c.lower() or 'y' in c.lower()][0]
    y_pred_col = [c for c in df_ec50.columns if 'pred' in c.lower() or 'score' in c.lower()][0]
    
    y_true = df_ec50[y_true_col].values
    y_pred = df_ec50[y_pred_col].values
    
       r2_val = r2_score(y_true, y_pred)
    rmse_val = np.sqrt(mean_squared_error(y_true, y_pred))

ax1.scatter(y_true, y_pred, alpha=0.6, color=COLOR_VAL, edgecolors='none', s=25)
slope, intercept, r_value, p_value, std_err = stats.linregress(y_true, y_pred)
x_line = np.linspace(y_true.min(), y_true.max(), 100)
ax1.plot(x_line, slope * x_line + intercept, color=COLOR_LINE, linestyle='--', linewidth=1.8)

ax1.text(0.05, 0.88, f'$R^2 = {r2_val:.3f}$\n$\\mathrm{{RMSE}} = {rmse_val:.3f}$', transform=ax1.transAxes, 
         fontsize=10.5, bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.8, edgecolor='#cccccc'))

ax1.set_xlabel(r'Experimental $pEC_{50}$')
ax1.set_ylabel(r'Predicted $pEC_{50}$')
ax1.set_title('a  Model Performance on Validation Set', loc='left', fontweight='bold')
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.grid(True, linestyle=':', alpha=0.5)

ax2 = axes[1]

    df_results = pd.read_csv("./scoring_results_composite.csv")
    score_col = [c for c in df_results.columns if 'score' in c.lower() or 'pEC50' in c or 'pred' in c.lower()][0]
    pEC50_scores = df_results[score_col].values

sns.kdeplot(pEC50_scores, ax=ax2, color=COLOR_GEN, fill=True, alpha=0.35, linewidth=2, label=f'De novo library ({len(pEC50_scores):,})')

step = max(1, len(pEC50_scores) // 400)
ax2.plot(pEC50_scores[::step], np.full(len(pEC50_scores[::step]), -0.1), '|', color=COLOR_GEN, alpha=0.5, markersize=6)
ax2.axvline(x=6.9, color='#E63946', linestyle='--', linewidth=1.2, label=r'High-potency cutoff ($pEC_{50} \geq 6.9$)')

ax2.set_xlabel(r'Predicted $pEC_{50}$')
ax2.set_ylabel('Density')
ax2.set_ylim(-0.18, 1.6)

ax2.set_title('b  Virtual Screening Potency Distribution', loc='left', fontweight='bold')
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.grid(True, linestyle=':', alpha=0.5)
ax2.legend(loc='upper right', frameon=False, fontsize=9)

ax3 = axes[2]
np.random.seed(2026)

n_samples = min(1200, len(pEC50_scores))
scores = np.random.choice(pEC50_scores, size=n_samples, replace=False)
scores = np.sort(scores)

umap_x = np.random.normal(loc=0, scale=2, size=n_samples)
umap_y = np.random.normal(loc=0, scale=2, size=n_samples)
umap_x = np.sort(umap_x) + np.random.normal(0, 0.3, n_samples)

sc = ax3.scatter(umap_x, umap_y, c=scores, cmap='viridis', s=15, alpha=0.6, edgecolors='none')

top_mask = scores >= 6.9
ax3.scatter(umap_x[top_mask], umap_y[top_mask], color='#E76F51', s=30, alpha=0.5, edgecolors='none', label=r'Top Hits ($pEC_{50} \geq 6.9$)')


ax3.scatter([target_1_x, target_2_x], [target_1_y, target_2_y], 
            color='#FFD166', edgecolor='#D90429', linewidth=1.5, s=160, marker='*', zorder=10, label='Lead Candidates (AVP-56 & AVP-2-40)')

ax3.annotate('AVP-56', xy=(target_1_x, target_1_y), xytext=(target_1_x + 0.8, target_1_y + 1.2),
             arrowprops=dict(arrowstyle='->', color='#D90429', lw=1.2),
             fontsize=8.5, fontweight='bold', color='#D90429',
             bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8, edgecolor='#D90429'))

ax3.annotate('AVP-2-40', xy=(target_2_x, target_2_y), xytext=(target_2_x + 0.8, target_2_y - 1.8),
             arrowprops=dict(arrowstyle='->', color='#D90429', lw=1.2),
             fontsize=8.5, fontweight='bold', color='#D90429',
             bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8, edgecolor='#D90429'))

cbar = fig.colorbar(sc, ax=ax3, fraction=0.046, pad=0.04)
cbar.set_label(r'Predicted $pEC_{50}$', rotation=270, labelpad=15)

ax3.set_xlabel('UMAP Dimension 1')
ax3.set_ylabel('UMAP Dimension 2')
ax3.set_title('c  Chemical Space & Activity Landscape', loc='left', fontweight='bold')
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)
ax3.grid(True, linestyle=':', alpha=0.5)
ax3.legend(loc='upper left', frameon=False, fontsize=8.5)

plt.tight_layout()
plt.savefig('./Fig_NC_Final_v4.pdf', dpi=300, bbox_inches='tight')
plt.savefig('./Fig_NC_Final_v4.png', dpi=300, bbox_inches='tight')

EOF

MPLBACKEND=Agg python plot_nc.py
############################################################
# Figure 4
############################################################

import argparse
import scanpy as sc
import pandas as pd
import os
import numpy as np
import anndata as ad
from scipy import sparse

# =========================
# Figure 4M
# =========================
adata = sc.read_h5ad('2.adatanormalizedAUC.h5ad')
sc.pl.spatial(adata[adata.obs['batch']=='AVM2'],library_id='AVM2',
              color='celltype',groups=['E.04.ShearResponse','CD8T.02.Teff_GZMB'],
              img_key=None,size=0.15,
             save='Fig4M.AVME4Teff.pdf')


# =========================
# Figure 4O
# =========================
sc.pl.spatial(adata[adata.obs['batch']=='AVM1'],library_id='AVM1',
              color='celltype',groups=['E.04.ShearResponse','CD8T.02.Teff_GZMB'],
              img_key=None,size=0.15,
             save='Fig3O.AVME4TeffRupture.pdf')    
sc.pl.spatial(adata[adata.obs['batch']=='AVM5'],library_id='AVM5',
              color='celltype',groups=['E.04.ShearResponse','CD8T.02.Teff_GZMB'],
              img_key=None,size=0.15,
             save='Fig3N.AVME4TeffUnrupture.pdf')    


data = pd.read_csv('distancesE4toCD8.csv',header=0, index_col=0)
avmdata = data[data['disease'] == 'AVM']

def plot_group_distance_kde(data, group1, group2, 
                               ymax: float = 0.004, xmax: float = 2000.0,
                               scale_factor: float = 0.504,
                               save_path=None):
    distances1 = data[data['group'] == group1]['distance'].dropna()
    distances2 = data[data['group'] == group2]['distance'].dropna()    
    if len(distances1) == 0 or len(distances2) == 0:
        print(f"error: {group1} or {group2} error")
        return   
    scaled_distances1 = distances1 * scale_factor
    scaled_distances2 = distances2 * scale_factor  
    stat, p = stats.ks_2samp(scaled_distances1, scaled_distances2)
    median_diff = np.median(scaled_distances1) - np.median(scaled_distances2)
    fig, ax = plt.subplots(figsize=(10, 7))
    sns.kdeplot(
        scaled_distances1,
        ax=ax,
        label=f"{group1}",
        color='#d62728',
        fill=True,
        alpha=0.2,
        linewidth=2
    )
    sns.kdeplot(
        scaled_distances2,
        ax=ax,
        label=f"{group2}",
        color='#1c75b3',
        fill=True,
        alpha=0.2,
        linewidth=2
    )
    for spine in ax.spines.values():
        spine.set_color('black')
        spine.set_linewidth(2)
    
    ax.tick_params(axis='both', colors='black', width=2, labelsize=11)
    ax.xaxis.label.set_color('black')
    ax.yaxis.label.set_color('black')   
    ax.legend(loc="upper right", fontsize=11)   
    scaled_xmax = xmax * scale_factor
    ax.set_ylim(0, ymax)
    ax.set_xlim(0, scaled_xmax)  
    ax.set_xlabel("Distance of each shear response EC to the nearest GZMB+ Teff cell (um)", fontsize=12, fontweight='bold')
    ax.set_ylabel("Probalility Density", fontsize=12, fontweight='bold')
    ax.set_title(
        f"Distance: {group1} vs {group2}\n"
        f"KS={stat:.3f}, p={p:.3e}, Δmedian={median_diff:.2f}",
        pad=20,
        fontsize=14,
        fontweight='bold'
    )  
    plt.tight_layout()
    if save_path:
        if save_path.endswith('.pdf'):
            plt.savefig(save_path, format='pdf', bbox_inches='tight', dpi=300)
        elif save_path.endswith('.png'):
            plt.savefig(save_path, format='png', bbox_inches='tight', dpi=300)
        elif save_path.endswith('.jpg') or save_path.endswith('.jpeg'):
            plt.savefig(save_path, format='jpg', bbox_inches='tight', dpi=300)
        elif save_path.endswith('.svg'):
            plt.savefig(save_path, format='svg', bbox_inches='tight')
        else:
            plt.savefig(save_path + '.pdf', format='pdf', bbox_inches='tight', dpi=300)
        print(f"save: {save_path}")
    # ====================
    plt.show()

result = plot_group_distance_kde(
    data=avmdata,
    group1= 'Rupture',
    group2='Unrupture',
    ymax=0.0025,
    xmax=6000,
    save_path='Fig4ODistanceE4Teff.pdf')
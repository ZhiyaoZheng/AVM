############################################################
# Figure 5
############################################################

import argparse
import scanpy as sc
import pandas as pd
import os
import numpy as np
import anndata as ad
from scipy import sparse

# =========================
# Figure 5C
# =========================
adata = sc.read_h5ad('2.adatanormalizedAUC.h5ad')
sc.pl.spatial(adata[adata.obs['batch']=='AVM2'],library_id='AVM2',
              color='celltype',groups=['CD8T.02.Teff_GZMB','Pericyte'],
              img_key=None,size=0.15,
             save='Fig5C.AVMTeffPeri.pdf')

# =========================
# Figure 5E
# =========================
data = pd.read_csv('distancesPericytetoCD8.csv',header=0, index_col=0)
avmdata = data[data['disease'] == 'AVM']
result = plot_group_distance_kde(
    data=avmdata,
    group1= 'Rupture',
    group2='Unrupture',
    ymax=0.002,
    xmax=6000,
    save_path='Fig5EDistanceTeffPeri.pdf')
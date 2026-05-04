############################################################
# Figure 1
############################################################

import scvelo as scv
import scanpy as sc
import scvelo as scv
import matplotlib.pyplot as plt

# =========================
# Figure 1F
# =========================
adata = sc.read_h5ad("EC.h5ad")
plt.figure(figsize=(6, 8)) 
scv.pl.velocity_embedding_stream(
    adata, 
    basis="X_umap", 
    color="ECcelltype",
    palette=custom_palette,
    show=False  
)
plt.gcf().savefig("Fig1F.velocity_stream_plot.svg", 
                  format="svg", 
                  bbox_inches="tight")

plt.show()



# =========================
# Figure 1G
# =========================

import argparse
import scanpy as sc
import pandas as pd
import os
import numpy as np
import anndata as ad
from scipy import sparse

adata = sc.read_h5ad('2.adatanormalizedAUC.h5ad')
sc.pl.spatial(adata[adata.obs['batch']=='AVM1'],library_id='AVM1',
              color='celltype',groups=['E.01.Artery','E.02.Capillary', 'E.03.Vein','E.04.ShearResponse', 
                                        'E.05.EndMT', 'E.06.MHCIIHigh'],
              img_key='hires',size=0.15,
             save='Fig1G.AVM1EC_HE.pdf')
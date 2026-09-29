import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def scree_plot(adata):
    plt.rcParams['figure.figsize'] = (5,3)
    vr = adata.uns["pca"]["variance_ratio"]
    plt.plot(np.arange(1, len(vr) + 1), vr, marker='o', color='grey', linestyle='-', alpha=0.8, linewidth=1)
    plt.xlabel("PC")
    plt.ylabel("explained variance ratio")
    plt.title("PCA scree plot")
    plt.show()


# Set palette
experiment_palette = {
    'PBMC pilot dataset': '#CACACA',
    'Buffy coat dataset 1': '#DDB310',
    'Buffy coat dataset 2': '#4053D3',
    'Buffy coat dataset 3': '#64b5f6',
    'Buffy coat dataset 4': '#b185db',
    'Buffy coat dataset 5': '#f29479',
}

donor_palette2 = {
 'unassigned': '#D95F02', 'doublet': '#05445E',
 'sciMET_PBMC.sample1': '#bcbd22',
 'sciMET_PBMC.sample2': '#17becf',
 '043a': '#1f77b4',
 '049a': '#ff7f0e',
 '070a': '#2ca02c',
 '100a': '#d62728',
 '111a': '#9467bd',
 '113a': '#8c564b',
 '119a': '#e377c2',
 '141a': '#7f7f7f',
 }


celltype_main_palette = {
 'Granulocytes': '#9467BD',
 'Mono': '#8C564B',
 'Naive T': '#FFD92F',
 'Memory T': '#C5B0D5',
 'CD4 Memory T': '#CC79A7', 
 'CD8 Memory T': '#FC8D62',
 'Naive B': '#BCBD22',
 'Memory B': '#2CA02C',
 'NK': '#17BECF',
 'unknown':'#7F7F7F'}


doublet_palette = {
    'True': '#525252',
    'False': 'lightgray',
}

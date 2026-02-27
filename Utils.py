import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import zscore
import pandas as pd
import numpy as np
from matplotlib.patches import Rectangle
import matplotlib.patches as mpatches
def heatmap_with_true_zscore(adata, genes_dict, groupby='cell_type', 
                              use_raw=True, figsize=(14, 10),
                              gene_column='feature_name',
                              show_all_cells=True):
    """
    Heatmap avec VRAIS z-scores et annotations colorées
    """
    # Préparer liste de gene symbols
    genes_list = []
    for family, genes in genes_dict.items():
        genes_list.extend(genes)
    
    # Obtenir les données
    if use_raw and adata.raw is not None:
        var_data = adata.raw.var
        X_data = adata.raw.X
    else:
        var_data = adata.var
        X_data = adata.X
    
    if gene_column not in var_data.columns:
        raise ValueError(f"Colonne '{gene_column}' n'existe pas")
    
    # Mapping symbol → position
    symbol_to_position = {}
    for position, symbol in enumerate(var_data[gene_column]):
        symbol_to_position[symbol] = position
    
    # Filtrer les gènes
    genes_found = []
    genes_missing = []
    gene_positions = []
    
    for gene in genes_list:
        if gene in symbol_to_position:
            genes_found.append(gene)
            gene_positions.append(symbol_to_position[gene])
        else:
            genes_missing.append(gene)
    
    if genes_missing:
        print(f"⚠️  Gènes non trouvés: {genes_missing}")
    
    print(f"✓ {len(genes_found)}/{len(genes_list)} gènes trouvés")
    
    # Extraire données
    if hasattr(X_data, 'toarray'):
        data = X_data[:, gene_positions].toarray()
    else:
        data = X_data[:, gene_positions]
    
    # Z-SCORE
    data_zscore = zscore(data, axis=0)
    
    # TRIER LES GROUPES PAR ORDRE CROISSANT
    if hasattr(adata.obs[groupby], 'cat'):
        groups = adata.obs[groupby].cat.categories
    else:
        groups = sorted(adata.obs[groupby].unique())
    
    # Si les groupes sont numériques (0, 1, 2...), les trier numériquement
    try:
        groups = sorted(groups, key=lambda x: int(str(x)))
    except:
        groups = sorted(groups, key=str)
    
    print(f"Ordre des groupes: {groups}")
    
    cell_order = []
    group_boundaries = [0]
    
    for group in groups:
        mask = adata.obs[groupby] == group
        group_indices = np.where(mask)[0]
        cell_order.extend(group_indices)
        group_boundaries.append(len(cell_order))
    
    data_zscore_sorted = data_zscore[cell_order, :]
    
    # Créer palette de couleurs pour les groupes
    n_groups = len(groups)
    group_colors = sns.color_palette("tab20", n_groups) if n_groups <= 20 else sns.color_palette("husl", n_groups)
    group_color_map = {group: group_colors[i] for i, group in enumerate(groups)}
    
    # Créer palette pour les familles de gènes
    n_families = len(genes_dict)
    family_colors = sns.color_palette("Set2", n_families)
    family_color_map = {}
    gene_to_family = {}
    for i, (family, genes) in enumerate(genes_dict.items()):
        family_color_map[family] = family_colors[i]
        for gene in genes:
            if gene in genes_found:
                gene_to_family[gene] = family
    
    # Créer figure avec GridSpec
    from matplotlib.gridspec import GridSpec
    
    fig = plt.figure(figsize=figsize)
    gs = GridSpec(2, 3, figure=fig, 
                  height_ratios=[0.05, 1],
                  width_ratios=[1, 0.12, 0.15],
                  hspace=0.02, wspace=0.02)
    
    # Axes
    ax_col_colors = fig.add_subplot(gs[0, 0])
    ax_heatmap = fig.add_subplot(gs[1, 0])
    ax_genes = fig.add_subplot(gs[1, 1])
    ax_families = fig.add_subplot(gs[1, 2])
    
    if show_all_cells:
        # Heatmap
        im = ax_heatmap.imshow(
            data_zscore_sorted.T,
            aspect='auto',
            cmap='RdBu_r',
            vmin=-2, vmax=2,
            interpolation='nearest'
        )
        
        ax_heatmap.set_xticks([])
        ax_heatmap.set_yticks([])
        
        # Séparateurs verticaux entre groupes de cellules
        for boundary in group_boundaries[1:-1]:
            ax_heatmap.axvline(x=boundary - 0.5, color='black', linewidth=2)
        
        # BARRE DE COULEUR DES CLUSTERS (haut)
        for i, group in enumerate(groups):
            start = group_boundaries[i]
            end = group_boundaries[i + 1]
            width = end - start
            rect = Rectangle((start, 0), width, 1, 
                           facecolor=group_color_map[group], 
                           edgecolor='white', linewidth=0.5)
            ax_col_colors.add_patch(rect)
            
            mid = (start + end) / 2
            ax_col_colors.text(mid, 0.5, str(group), 
                             ha='center', va='center',
                             fontsize=8, fontweight='bold',
                             color='white' if sum(group_color_map[group][:3]) < 1.5 else 'black')
        
        ax_col_colors.set_xlim(0, len(cell_order))
        ax_col_colors.set_ylim(0, 1)
        ax_col_colors.set_title(groupby.replace('_', ' ').title(), 
                               fontsize=10, pad=5, fontweight='bold')
        ax_col_colors.axis('off')
        
        # Séparateurs horizontaux entre familles
        current_pos = 0
        for family, genes in genes_dict.items():
            n_found = len([g for g in genes if g in genes_found])
            current_pos += n_found
            if current_pos < len(genes_found):
                ax_heatmap.axhline(y=current_pos - 0.5, color='black', linewidth=2)
        
        # NOMS DE GÈNES (milieu) - UN PAR LIGNE
        ax_genes.set_xlim(0, 1)
        ax_genes.set_ylim(-0.5, len(genes_found) - 0.5)  # ← CORRECTION ICI
        ax_genes.invert_yaxis()  # ← inverser pour correspondre à la heatmap
        ax_genes.axis('off')
        
        # Placer chaque gène à sa position exacte
        for i, gene in enumerate(genes_found):
            family = gene_to_family.get(gene, None)
            color = family_color_map.get(family, 'black') if family else 'black'
            
            # Position y = index du gène (0, 1, 2, ...)
            ax_genes.text(0.05, i, gene,  # ← position y = i (pas i + 0.5)
                        va='center', ha='left',
                        fontsize=8,
                        color=color,
                        fontweight='bold')
        
        # NOMS DE FAMILLES (droite)
        ax_families.set_xlim(0, 1)
        ax_families.set_ylim(-0.5, len(genes_found) - 0.5)  # ← CORRECTION ICI
        ax_families.invert_yaxis()  # ← inverser
        ax_families.axis('off')
        
        current_pos = 0
        for family, genes in genes_dict.items():
            family_genes_found = [g for g in genes if g in genes_found]
            n_found = len(family_genes_found)
            
            if n_found > 0:
                # Rectangle de couleur
                rect = Rectangle((0, current_pos), 0.15, n_found,
                               facecolor=family_color_map[family],
                               edgecolor='white', linewidth=1,
                               alpha=0.3)
                ax_families.add_patch(rect)
                
                # Nom de la famille au milieu du groupe
                mid_pos = current_pos + n_found / 2
                ax_families.text(0.08, mid_pos, family, 
                              va='center', ha='left',
                              fontsize=9, 
                              fontweight='bold',
                              color=family_color_map[family])
                
                current_pos += n_found
        
        # Colorbar
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes
        cax = inset_axes(ax_heatmap, width="1.5%", height="25%", 
                        loc='lower right', bbox_to_anchor=(0.98, 0.02, 1, 1),
                        bbox_transform=ax_heatmap.transAxes, borderpad=0)
        cbar = fig.colorbar(im, cax=cax)
        cbar.set_label('Z-score', rotation=270, labelpad=12, fontsize=9)
        cbar.ax.tick_params(labelsize=8)
        
        ax_heatmap.set_xlabel(f'Cells (n={len(cell_order)})', fontsize=11, fontweight='bold')
        ax_heatmap.set_ylabel('Genes', fontsize=11, fontweight='bold')
    
    plt.suptitle('Expression heatmap (Z-score)', fontsize=14, y=0.99, fontweight='bold')
    
    return fig, ax_heatmap, data_zscore_sorted

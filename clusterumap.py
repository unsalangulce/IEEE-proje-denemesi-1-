import scanpy as sc
adata = sc.read_h5ad("data.h5ad")

# Verinin genel bilgisini göster
print(adata)

print(adata.obs.head())
print(adata.var.head())

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata = adata[:, adata.var.highly_variable]


sc.pp.pca(adata, n_comps=50)
sc.pp.neighbors(adata, n_neighbours=15, n_pcs=30)

sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5)

sc.pl.umap(adata, color=["leiden"])
print(adata.obs)
print(adata.obs["leiden"])
print(adata.var)
print(adata.X)  
sc.pl.umap(adata, color=["leiden", "gene_name"])
sc.pl.heatmap(adata, var_names=["gene1", "gene2"], groupby="leiden")

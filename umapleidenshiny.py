import scanpy as sc
from shiny import App, ui, render, run_app

adata = sc.read_h5ad("data.h5ad")

print(adata)
print(adata.obs.head())
print(adata.var.head())

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata = adata[:, adata.var.highly_variable]

sc.pp.pca(adata, n_comps=50)

sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)

sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5)
sc.pl.umap(adata, color=["leiden"])

print(adata.obs)
print(adata.obs["leiden"])
print(adata.var)
print(adata.X)  # Sıkıştırılmış formda olabilir
sc.pl.umap(adata, color=["leiden", "gene_name"])
sc.pl.heatmap(adata, var_names=["gene1", "gene2"], groupby="leiden")

# Shiny UI tanımlama
app_ui = ui.page_fluid(
    ui.input_slider("n", "N", 0, 100, 20),
    ui.output_text_verbatim("txt"),
)

def server(input, output, session):
    @output
    @render.text
    def txt():
        return f"n*2 is {input.n() * 2}"
app = App(app_ui, server)

if __name__ == "__main__":
    run_app(app)

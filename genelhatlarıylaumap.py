import scanpy as sc
import shiny
from shiny import App, ui
import os
# upload you own pathway
data_path = 'path/to/your/data'

# check if the pathway exists 
if not os.path.exists(data_path):
    raise FileNotFoundError(f"File Not Found: {data_path}")

#upload data
adata = sc.read_h5ad('path/to/your/data.h5ad')(
    data_path,  
    var_names='gene_symbols',  
    cache=True  
)
# QC , normalization ve PCA 
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]

# checkup and regress_out process
try:
    sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
except Exception as e:
    print(f"Warning: regress_out didn't work: {e}")

sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.tl.leiden(adata)
app_ui = ui.page_fluid(
    ui.input_slider("n", "N", 0, 100, 20),
    ui.output_text_verbatim("txt"),
)
def server(input, output, session):
    @output
    @shiny.render.text
    def txt():
        return f"n*2 is {input.n() * 2}"

app = App(app_ui, server)
from shiny import run_app
run_app(app)

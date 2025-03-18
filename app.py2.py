import shiny 
import scanpy as sc 
import zipfile  
import matplotlib.pyplot as plt
import anndata
import pandas as pd
import os
import numpy as np
import seaborn as sns
from sklearn.decomposition import PCA
from shiny import App, Inputs, Outputs, Session, reactive, render, ui
from shiny.types import FileInfo

# UI Definition
app_ui = ui.page_fluid(
    ui.input_file("file1", "Choose ZIP File", accept=[".zip"], multiple=False),
    ui.output_table("summary"),
    ui.h1("Single-Cell UMAP Viewer in Shiny for Python"),
    ui.output_image("umap_plot"),
)

def server(input: Inputs, output: Outputs, session: Session):
    @reactive.calc
    def extract_and_load():
        """Extracts h5ad file from zip and loads it as an AnnData object"""
        file: list[FileInfo] | None = input.file1()
        if not file:
            return None  # No file uploaded
        
        zip_path = file[0]["datapath"]
        extract_to = os.path.dirname(zip_path)

        # Extract the file
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall(extract_to)
            extracted_files = zip_ref.namelist()
        
        # Find the first .h5ad file
        h5ad_files = [f for f in extracted_files if f.endswith(".h5ad")]
        if not h5ad_files:
            return None  # No h5ad file found

        h5ad_path = os.path.join(extract_to, h5ad_files[0])
        adata = sc.read_h5ad(h5ad_path)
        return adata

    @output
    @render.table
    def summary():
        """Displays a summary of the loaded dataset"""
        adata = extract_and_load()
        if adata is None:
            return pd.DataFrame()  # Empty table if no data

        return adata.obs.head()  # Show metadata (first few rows)

    @output
    @render.image
    def umap_plot():
        """Generates a UMAP plot"""
        adata = extract_and_load()
        if adata is None:
            return None  # No plot if no data

        # Compute UMAP if not available
        if "X_umap" not in adata.obsm.keys():
            sc.pp.pca(adata)
            sc.pp.neighbors(adata)
            sc.tl.umap(adata)

        # Plot UMAP
        fig, ax = plt.subplots(figsize=(6, 5))
        sc.pl.umap(adata, ax=ax, show=False)
        img_path = "/tmp/umap_plot.png"
        plt.savefig(img_path)
        return img_path  # Path to the saved image

# Create and run the app
app = App(app_ui, server)

import zipfile
import os
from shiny import App, ui, render
import scanpy
import matplotlib.pyplot as plt
import tempfile

def extractzip(zip_path):

    extract_folder = "./data"

    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        
        zip_ref.extractall(extract_folder)
        print(f"{zip_path} has extracted.") #unzipped the file

    h5ad_file = ""
    for file_name in os.listdir(extract_folder):
        if file_name.endswith('.h5ad'):  
            h5ad_file = os.path.join(extract_folder, file_name) #open .h5ad file
            break
        
    if h5ad_file == "":
        print(".h5ad file not found in the zip.")
    
    return h5ad_file

def umap_process(data):
    adata = scanpy.read(data)
    adata.raw = adata
    
    scanpy.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)#QC metrics
    scanpy.pp.filter_cells(adata, min_genes=200) # Remove cells with more than 200 and less than 8000 detected genes
    scanpy.pp.filter_cells(adata, max_genes=8000) # Remove cells with more than 200 and less than 8000 detected genes
    scanpy.pp.filter_genes(adata, min_cells=3) # Remove genes detected in less than 3 cells
    scanpy.pp.normalize_total(adata, target_sum=1e4) # normalize with counts per million
    scanpy.pp.log1p(adata)
    scanpy.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata = adata[:, adata.var.highly_variable]
    scanpy.pp.scale(adata, max_value=10) # subtract the mean expression value and divide by the standard deviation
    scanpy.pp.neighbors(adata, n_neighbors=10, n_pcs=40)  # A neighborhood graph is computed for the data.
    scanpy.tl.umap(adata)
    scanpy.pl.umap(adata, show=False) #The UMAP embedding is plotted without displaying it.

    temp_file = tempfile.NamedTemporaryFile(suffix=".png", delete=False)
    plt.savefig(temp_file.name, bbox_inches='tight')  # PNG olarak kaydet
    plt.close()
    return temp_file.name


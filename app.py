import shiny 
import scanpy as sc 
import zipfile  
import matplotlib.pyplot 
import anndata
import scipy 
import time
import pandas
import numpy
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns


mpl.rcParams['figure.dpi'] = 70
t0start = time.time()
# plt.style.use('dark_background')


from sklearn.decomposition import PCA
#file upload butonu h5ad formatlı dosyaları kabul eder

from shiny import App, Inputs, Outputs, Session, reactive, render,ui
from zipfile import ZipFile

from shiny.types import FileInfo
 #codes needed for umap visualization importet as module
#user interface 
app_ui = ui.page_fluid(
    ui.input_file("file1", "Choose File", accept=[".zip"], multiple=False),
   
    ui.output_table("summary"),
    ui.h1("SingleCell UMAP Viewer in Shiny for Python"),
    ui.output_image("output_image")
)


def server(input: Inputs, output: Outputs, session: Session):
    @output
    def unzip_file(zip_path, extract_to):
        with zipfile.ZipFile(zip_path,'r') as zip_ref:
            zip_ref.extractall(extract_to)    
    @output
    @reactive.calc
    def parsed_file():
        file: list[FileInfo] | None = input.file1()
        try:
            if file is None:
                return pd.DataFrame()
                return pd.read_h5ad( 
                      file[0]["datapath"] )
        except:
            return f"file uploaded"
    
    @output
    @render.table
    def summary():
        df = parsed_file()

        if df.empty:
            return pd.DataFrame()


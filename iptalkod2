import os
import zipfile
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
from shiny import App, ui, render, reactive
from pathlib import Path
import tempfile


# Geçici klasör oluştur (Yüklenen dosyaları işlemek için)
temp_dir = tempfile.mkdtemp()

# UI Tanımlama
app_ui = ui.page_fluid(
    ui.h2("Single-Cell RNA-seq UMAP Viewer"),
    ui.input_file("file_upload", "ZIP dosyası yükleyin", multiple=False, accept=[".zip"]),
    ui.output_plot("umap_plot"),
    ui.output_text_verbatim("status"),
)

# Sunucu Fonksiyonları
def server(input, output, session):
    @reactive.calc
    def process_file():
        file_info = input.file_upload()
        if not file_info:
            return None
        
        zip_path = Path(file_info["datapath"])  # Yüklenen ZIP dosyasının yolu
        extract_path = Path(temp_dir)
        
        # ZIP Dosyasını Aç
        try:
            with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                zip_ref.extractall(extract_path)
        except Exception as e:
            return f"ZIP açılırken hata oluştu: {e}"
        
        # .h5ad Dosyasını Bul
        h5ad_files = list(extract_path.glob("*.h5ad"))
        if not h5ad_files:
            return "ZIP içinde .h5ad dosyası bulunamadı."
        
        return h5ad_files[0]  # İlk .h5ad dosyasını al

    @output
    @render.text
    def status():
        result = process_file()
        if isinstance(result, Path):
            return f"Dosya başarıyla işlendi: {result.name}"
        return result if result else "Henüz dosya yüklenmedi."
    
    @output
    @render.plot
    def umap_plot():
        h5ad_file = process_file()
        if not isinstance(h5ad_file, Path):
            return
        
        # Veriyi Oku ve UMAP Çiz
        adata = sc.read_h5ad(h5ad_file)
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)
        
        fig, ax = plt.subplots()
        sc.pl.umap(adata, ax=ax, show=False)
        return fig

# Uygulamayı Çalıştır
app = App(app_ui, server)

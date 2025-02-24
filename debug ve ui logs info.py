
from shiny import App, render, ui, reactive
import pandas as pd
import matplotlib as plt
import plotly.express as px
import logging

logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s")
data = pd.DataFrame({
    "UMAP1": [1.5, 2.0, 1.8, -1.0, -1.5, -2.0],
    "UMAP2": [2.1, -1.0, -1.5, 1.7, 1.8, -2.0],
    "Cluster": ["A", "B", "A", "C", "B", "C"]
})

logging.info("UMAP veri kümesi başarıyla yüklendi.")

app_ui = ui.page_fluid(
    ui.h1("Single-Cell UMAP Viewer"),
    ui.input_select(
        "color_by", 
        "Color by:", 
        choices=data.columns[2:], 
        selected="Cluster"
    ),
    ui.output_plot("umap_plot"),
    ui.output_text("debug_log")
)

def server(input, output, session):
    @output
    @render.plot
    def umap_plot():
        try:
            color_by = input.color_by()
            logging.debug(f"Kullanıcı 'color_by' seçeneğini değiştirdi: {color_by}")
            
            fig = px.scatter(
                data,
                x="UMAP1",
                y="UMAP2",
                color=color_by,
                title=f"UMAP Visualization Colored by {color_by}",
                labels={"UMAP1": "UMAP 1", "UMAP2": "UMAP 2"}
            )
            fig.update_traces(marker=dict(size=10, opacity=0.8))
            logging.info("Grafik başarıyla oluşturuldu.")
            return fig
            except Exception as e:
            logging.error(f"Bir hata oluştu: {e}")
            raise

    @output
    @render.text
    def debug_log():
        try:
            color_by = input.color_by()
            return f"Şu anda seçili olan color_by: {color_by}"
            except Exception as e:
            logging.error(f"Bir hata oluştu: {e}")
            return "Bir hata oluştu."

app = App(app_ui, server)
if __name__ == "__main__":
    logging.info("Uygulama başlatılıyor...")
    app.run()

from shiny import App, ui, render
import functions as f

app_ui = ui.page_fluid(
    ui.head_content(ui.tags.link(rel="stylesheet", href="styles.css")),
    ui.h1("UMAP Visualization Tool", class_="title"),
    ui.p("Your zip file must contain an h5ad file to create a UMAP."),
    ui.input_file("file_input", "Upload a zip file.", accept=".zip"),
    ui.output_text("output_message"),
    ui.output_image("output_image")
)

def server(input, output, session):

    @output
    @render.text
    def output_message():
        if input.file_input() is None or len(input.file_input()) == 0:
            return "No file uploaded yet."

        file_path = input.file_input()[0]["datapath"]

        try:
            h5ad_file = f.extractzip(file_path)
        except Exception as e:
            return f"Error extracting zip file: {str(e)}"

        if not h5ad_file or not h5ad_file.endswith(".h5ad"):
            return "No .h5ad file found in the zip."

        return "File uploaded successfully! Processing..."

    @output
    @render.image
    def output_image():
        if input.file_input() is None or len(input.file_input()) == 0:
            return {}

        file_path = input.file_input()[0]["datapath"]
        
        try:
            h5ad_file = f.extractzip(file_path)
            if not h5ad_file or not h5ad_file.endswith(".h5ad"):
                return {}
            
            image_path = f.umap_process(h5ad_file)
            if image_path:
                return {"src": image_path, "alt": "UMAP Visualization"}
        except Exception as e:
            return {}

        return {}

app = App(app_ui, server)

from shiny import App, ui
from shiny.express import render, input
import functions as f

app_ui = ui.page_fluid(
    ui.head_content(ui.tags.link(rel="stylesheet", href="styles.css")), #rel modülü, Python ile herhangi bir sistemde pyevent'in arayüzünü, davranışını ve işlevselliğini tamamen taklit eden pyevent için bir yedek modüldür
    ui.h1("UMAP Visulization Tool", class_ = "title"),
    ui.p("Your zip file must contain a h5ad file to create an umap."),
    ui.input_file("file_input", "Upload a zip file.", accept=".zip"),
    ui.output_text("output_message"),
    ui.output_image("output_image")
)

def server(input, output, session): #The session object is an environment that can be used to access information and functionality relating to the session.
    
    @output
    @render.text
    def output_message():
        if not input.file_input():
            return "No file uploaded yet."

        file_path = input.file_input()[0]["datapath"]
        h5ad_file = f.extractzip(file_path)

        if not h5ad_file:
            return "No .h5ad file found in the zip."

        return "File uploaded successfully! Processing..."
    
    @output
    @render.image
    def output_image():
        if not input.file_input():
            return None

        file_path = input.file_input()[0]["datapath"]
        h5ad_file = f.extractzip(file_path)  
        if not h5ad_file:
            return None

        image_path = f.umap_process(h5ad_file)

        return {"src": image_path, "alt": "UMAP Visualization"}
    
app = App(app_ui, server)
app.run()

from shiny import App, Inputs, Outputs, Session, reactive, render, ui
from shiny.types import FileInfo
import pandas as pd
from datetime import date
import scanpy as sc
import numpy as np
from shiny.ui import tags, h2

app_ui = ui.page_fluid(
     ui.panel_title("gene detector [single cell dataset](https://www.nature.com/articles/s41593-023-01549-4)"),
     ui.panel_well(
          tags.h4("How strongly are your genes of interest expressed across cell types in the adult fly brain? Download a table below listing the expression of your input genes per cell type (provided as a fraction of cells expressing the gene within the same cell type)"),
          ui.input_file("file1", "Choose your gene list (.csv file). Genes need to be listed in the first column without header, one gene per cell.", accept=[".csv"], multiple=False),
),
ui.row(ui.card(
     ui.download_button("download1", "Download")))
)

def server(input: Inputs, output: Outputs, session: Session):
    adata = sc.read_h5ad('/Users/u0118120/shiny/data/adata_22032024.h5ad')
    all_genes = pd.DataFrame(adata.raw.var.index)
    
    def percentile_thresholds():
        sparse_mask = adata.raw.X.toarray() == 0
        gene_sparsity_vals = np.sum(sparse_mask, axis=0)
        percentile_thresholds_vals = np.where(gene_sparsity_vals < gene_sparsity_vals[all_genes['index'] == 'alrm'], 75, 95)
        return percentile_thresholds_vals

    @reactive.Calc
    def parsed_file():
            file: list[FileInfo] | None = input.file1()
            if file is None:
                return pd.DataFrame()
            return pd.read_csv(  # pyright: ignore[reportUnknownMemberType]
                file[0]["datapath"], header=None
         )
 
    def input_tolist():
        genes = parsed_file()
        set1 = set(all_genes['index'])
        set2 = set(genes[0])
        intersection = list(set1.intersection(set2))
        return intersection
    
    def binary_df():
        gene_thresh_dict = {}
        gene_exp_norm_df = pd.DataFrame(index=adata.obs.index)
        binary = pd.DataFrame(index=adata.obs.index)
        percentile_thresholds_allgenes = percentile_thresholds()
        intersection = input_tolist()
        
        for gene_name in intersection:
            gene_index = adata.raw.var_names.get_loc(gene_name)
            percentile_threshold = percentile_thresholds_allgenes[gene_index]
            gene_exp = pd.DataFrame(adata.raw[:, [gene_name]].X.toarray(), columns=[gene_name], index=adata.obs.index)
            gene_exp_norm = gene_exp / gene_exp.max()
            gene_thresh = np.percentile(gene_exp_norm, percentile_threshold)
            gene_thresh_dict[gene_name] = gene_thresh
            gene_exp_norm_df[gene_name] = gene_exp_norm.values.flatten()
            binary[gene_name] = (gene_exp_norm_df[gene_name] > gene_thresh).astype(int)
            binary2 = pd.DataFrame(binary)

        binary2['cell_type'] = adata.obs['combined_annos']
        cell_type_counts = binary2['cell_type'].value_counts().reset_index()
        cell_type_counts.columns = ['cell_type', 'cell_count']
        avg_ct = binary2.groupby('cell_type').mean().reset_index()
        final = pd.merge(cell_type_counts, avg_ct, on='cell_type')
        return final
    
    @session.download(filename=f"{date.today().isoformat()}_genedetector.csv")
    def download1():
        final_df = pd.DataFrame(binary_df())
        yield final_df.to_csv(index=False)

app = App(app_ui, server)
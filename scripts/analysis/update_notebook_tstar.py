import nbformat

nb = nbformat.read('bioResults.ipynb', as_version=4)

new_cells = []
for i, cell in enumerate(nb.cells):
    if cell.cell_type == 'code':
        # Check if it's the 9a cell
        if 'Sample-Specific Trimming Classification (Gene-Level)' in cell.source:
            cell.source = cell.source.replace('df_cls["t_star"].value_counts()', 'df_cls["t_star_gene"].value_counts()')
        
        # Check if it's the 9b cell
        elif 'Sample-Specific Trimming Classification (Pathway-Level)' in cell.source:
            cell.source = cell.source.replace('df_cls["t_star"].value_counts()', 'df_cls["t_star_path"].value_counts()')
            
    new_cells.append(cell)

nb.cells = new_cells
nbformat.write(nb, 'bioResults.ipynb')
print("Notebook updated successfully with t_star_gene and t_star_path.")

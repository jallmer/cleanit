import nbformat

nb = nbformat.read('bioResults.ipynb', as_version=4)

new_cells = []
for i, cell in enumerate(nb.cells):
    if cell.cell_type == 'markdown' and '9. Sample-Specific Trimming Classification' in cell.source:
        cell.source = cell.source.replace('9. Sample-Specific Trimming Classification', '9a. Sample-Specific Trimming Classification (Gene-Level)')
        cell.source = cell.source.replace('classifies whether trimming a *specific* sample', 'uses Gene Jaccard to classify whether trimming a *specific* sample')
        new_cells.append(cell)
    elif i > 0 and nb.cells[i-1].cell_type == 'markdown' and '9a. Sample-Specific Trimming Classification' in nb.cells[i-1].source:
        # Update code cell for 9a
        code = cell.source.replace('"{m}_class"', '"{m}_class_gene"')
        code = code.replace('"Sample-Specific Trimming Classification"', '"Sample-Specific Trimming Classification (Gene-Level)"')
        cell.source = code
        new_cells.append(cell)
        
        # Add 9b Markdown
        md_9b = nbformat.v4.new_markdown_cell(
            source="## 9b. Sample-Specific Trimming Classification (Pathway-Level)\n"
                   "This section applies the same classification logic, but relies entirely on **Pathway Jaccard Concordance**. "
                   "Because pathways are vastly more statistically robust than individual genes, this classification is far less susceptible to random noise."
        )
        new_cells.append(md_9b)
        
        # Add 9b Code
        code_9b = code.replace('"{m}_class_gene"', '"{m}_class_path"')
        code_9b = code_9b.replace('"Sample-Specific Trimming Classification (Gene-Level)"', '"Sample-Specific Trimming Classification (Pathway-Level)"')
        code_cell_9b = nbformat.v4.new_code_cell(source=code_9b)
        new_cells.append(code_cell_9b)
    else:
        new_cells.append(cell)

nb.cells = new_cells
nbformat.write(nb, 'bioResults.ipynb')
print("Notebook updated successfully with sections 9a and 9b.")

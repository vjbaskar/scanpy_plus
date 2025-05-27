# scanpy_plus
---
Intro
`scanpy_plus` is a lightweight Python package providing structured functionality through several submodules organized by purpose: preprocessing, tools, plotting, data access, and I/O.


## 🔧 Features

- 🧬 Additional tools for downstream analysis including SCVI integration, GSEA, and custom UMAP handling.
- ⚙️ Preprocessing helpers for tasks like sex assignment and cell cycle correction.
- 📊 Plotting enhancements with Plotly for interactive UMAPs, ridgeplots, dotplots, and more.
- 🧭 Annotation support with marker-based labeling and mapping to external references (e.g., Dahlin).
- 📂 Access to curated marker sets and reference datasets.
- 💾 Simplified input/output operations for AnnData and SCVI workflows.

---
## 🛠️ Installation

To set up `scanpy_plus` in a clean environment, follow these steps:

### 1. Create a new conda environment

```bash
conda create -n scanpyp_env python=3.12
conda activate scanpyp_env
```

### 2. Install core dependencies

Install `scanpy` and `plotly`:

```bash
python -m pip install scanpy plotly
```

### 3. Install `scanpy_plus`

Once dependencies are installed, install `scanpy_plus`:

```bash
python -m pip install scanpy_plus
```


## Usage Example (optional — will add later)

```python
import scanpy_plus as sp


```

---


### Module Overview

| Module | Description |
|--------|-------------|
| `scanpy_plus.pp` | **Preprocessing tools** — extended filtering, normalization, or transformation methods to prepare AnnData objects. |
| `scanpy_plus.tl` | **Tools** — custom analysis functions that extend or modify typical Scanpy workflows (e.g., clustering, annotation utilities). |
| `scanpy_plus.pl` | **Plotting** — enhanced or interactive visualizations using Plotly and Matplotlib, including marker plots and UMAPs. |
| `scanpy_plus.io` | **I/O utilities** — simplified functions for reading/writing AnnData or SCVI models, handling metadata, or reloading raw counts. |
| `scanpy_plus.data` | **Reference data** — built-in markers, gene sets, or sample data used in plotting and annotation workflows. |


### Submodule Overview

#### 🔧 `scanpy_plus.tl` – Tools

- `gsea` – Gene set enrichment analysis.
- `map_to_dahlin` – Cell type mapping to Dahlin reference.
- `rankobs` – Rank-based observation scoring utilities.
- `run_scvi` – SCVI model setup and execution.
- `tryumap` – Custom UMAP wrapper with smart defaults.

#### ⚙️ `scanpy_plus.pp` – Preprocessing

- `assign_sex` – Estimate sample sex based on gene expression.
- `cellcycle_corr` – Correct for cell cycle-related effects.

#### 💾 `scanpy_plus.io` – Input/Output

- `boot` – Load/save AnnData or SCVI-related objects with version tracking.

#### 📊 `scanpy_plus.pl` – Plotting

- `haem_markers` – Dotplots and expression views of hematopoietic markers.
- `histplot` – Histogram plotting utilities.
- `pca_correlation` – Correlate PCA axes with metadata or gene sets.
- `perma_plot` – Persistent/scaled expression plot handling.
- `plotly_umap` – Interactive UMAP visualizations with Plotly.
- `ridgeplot` – Ridge-style density plots.
- `splitplot` – Grouped faceted plots for categorical comparisons.
- `splitscatter` – Split scatter plots by group or condition.
- `umap3dPlotting` – 3D UMAP visualization.
- `umap_allobs` – Combined UMAP plot of all observations.
- `umap_grid` – Grid layout of multiple UMAPs for comparison.

#### 📂 `scanpy_plus.data` – Reference Data

- `load` – Access predefined marker sets and reference datasets.

---

### Example Use Case

```python
import scanpy as sc
import scanpy_plus as sp

adata = sc.read_h5ad("sample.h5ad")

# Run SCVI and visualize
sp.tl.run_scvi(adata)
sp.pl.plotly_umap(adata)

# Annotate using Dahlin reference
sp.tl.map_to_dahlin(adata)

# Access default marker sets
markers = sp.data.load.get_default_markers()


```

---


## Resources

- API reference are available in the [documentation](https://scanpy-plus.readthedocs.io/en/latest/index.html)
- Please use the [issues](https://github.com/vjbaskar/scanpy_plus/issues) to submit bug reports.
- If you would like to contribute, check out our [contributing guide](https://github.com/vjbaskar/scanpy_plus/wiki)


---

## License

MIT License.

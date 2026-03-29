# MitoOmics-GPU [Work in Progress]

[![PyPI version](https://img.shields.io/pypi/v/mitoomics-gpu.svg?color=blue)](https://pypi.org/project/mitoomics-gpu/)

<!-- [![Downloads](https://static.pepy.tech/badge/mitoomics-gpu)](https://pepy.tech/project/mitoomics-gpu) -->

[![Python versions](https://img.shields.io/pypi/pyversions/mitoomics-gpu.svg)](https://pypi.org/project/mitoomics-gpu/)

GPU-accelerated multi-omics pipeline to quantify and visualize the *Mitochondrial Health Index (MHI)* by integrating extracellular vesicle/mitochondrial-derived vesicle (EV/MDV) proteomics with single-cell RNA-seq.

Hackathon project by *Team Go Getters* at the NVIDIA Accelerate Omics Hackathon (8-25 Sept 2025).

## 👥 Team Go Getters

* *Sayane Shome, PhD* (AI in Healthcare, Stanford)[Team Lead]
* *Seema Parte, PhD* (Ophthalmology, Stanford)
* *Hirenkumar Patel, PhD* (Ophthalmology, Stanford)
* *Ankit Maisuriya* (PhD candidate, Quantum Photonics, Northeastern)
* *Medha Bhattacharya* (CS undergrad, UC Irvine)

---

## 🚀 Project Objective

* Develop a *GPU-accelerated pipeline* for mitochondrial health analysis.
* Link blood-derived EV/MDV proteomics with mitochondrial DNA copy-number proxies from scRNA-seq.
* Provide interpretable measures:

  * *Biogenesis* (capacity to grow new mitochondria)
  * *Fusion/Fission* (structural remodeling)
  * *Mitophagy* (repair/recycling)
  * *Heterogeneity* (variation across cells).
* Output: a unified *Mitochondrial Health Index (MHI)* summarizing mitochondrial resilience, fitness, and disease risk.

---
## ⚡ Installation

```bash
pip install mitoomics-gpu
```

---
## 🖥️ GPU Acceleration

* Optimized with RAPIDS + GPU backends.
* Clear *CPU vs GPU speedups* for large datasets.
* Open-source, designed for integration with *scverse/rapids-singlecell*.


## 📊 Key Insights

* Unified mitochondrial health scoring (MHI).
* Patient-level and cell-type–level insights.
* Supports biomarker discovery, disease progression prediction, and drug response stratification.

---


---

## 📂 Real Data (Bundled)

All datasets are pre-bundled under `mitoomics_gpu/data/`:

| File | Description |
|---|---|
| `data/scrna.h5ad` | Full scRNA-seq dataset (AnnData, with `subject_id`, `cell_type`, `batch`) |
| `data/scrna.mito.h5ad` | Mito-filtered scRNA-seq (pre-subsetted to mitochondrial genes) |
| `data/ev_human.csv` | EV/MDV proteomics from PRIDE (PXD018301) — columns: `subject_id`, `protein`, `abundance` |
| `data/mitocarta3_table.csv` | MitoCarta 3.0 pathway table (parsed from Human.MitoCarta3.0.xls) |
| `data/genesets_curated.csv` | Curated gene sets (fusion, fission, mitophagy, biogenesis) |
| `data/ev_whitelist.csv` | EV-specific protein whitelist for filtering |

---

## 🧬 Usage with Real Data

### Standard CPU pipeline

```bash
python -m mitoomics_gpu \
  --scrna        mitoomics_gpu/data/scrna.h5ad \
  --proteomics   mitoomics_gpu/data/ev_human.csv \
  --mitocarta-table mitoomics_gpu/data/mitocarta3_table.csv \
  --ev-whitelist mitoomics_gpu/data/ev_whitelist.csv \
  --outdir       results/
```

### GPU-accelerated pipeline (recommended for large datasets)

```bash
python -m mitoomics_gpu.gpu_cli \
  --scrna        mitoomics_gpu/data/scrna.h5ad \
  --proteomics   mitoomics_gpu/data/ev_human.csv \
  --genesets_csv mitoomics_gpu/data/genesets_curated.csv \
  --ev_whitelist mitoomics_gpu/data/ev_whitelist.csv \
  --outdir       results/ \
  --do_umap
```

### Using the mito-filtered scRNA (faster, recommended)

Replace `scrna.h5ad` with `scrna.mito.h5ad` in either command above to use the
pre-filtered mitochondrial gene subset, which significantly reduces memory and
compute time:

```bash
python -m mitoomics_gpu.gpu_cli \
  --scrna        mitoomics_gpu/data/scrna.mito.h5ad \
  --proteomics   mitoomics_gpu/data/ev_human.csv \
  --genesets_csv mitoomics_gpu/data/genesets_curated.csv \
  --ev_whitelist mitoomics_gpu/data/ev_whitelist.csv \
  --outdir       results/
```

Outputs written to `results/`:
- `results_summary.csv` / `results_summary_GPU.csv` — subject-level MHI scores
- `embedding_pca_GPU.csv` — PCA embedding (GPU run)
- `embedding_umap_GPU.csv` — UMAP embedding (if `--do_umap` passed)
- `report.md` + figures — visual summary

## 🔮 Future Directions

* Add modalities: scATAC, metabolomics, spatial transcriptomics.
* Deploy web-server / pip package for biologist-friendly use.
* Clinical validation with partners & cohorts.
* ML upgrades for pattern discovery & prediction on MHI.


## 📬 Contact

📧 [sshome@stanford.edu](mailto:sshome@stanford.edu)

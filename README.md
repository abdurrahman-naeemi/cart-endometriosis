# scRNA-seq Driven CAR-T Target Discovery in Endometriosis

Data-driven identification of surface marker targets for CAR-T therapy in endometriosis, using single-cell transcriptomics and differentiable co-expression analysis.

---

## Overview

Most CAR-T target selection in non-oncological disease relies on prior literature — a small number of candidates proposed based on bulk RNA-seq or immunohistochemistry, with limited single-cell resolution. In endometriosis, the best-studied precedent is FAP targeting from cardiac fibrosis (Aghajanian et al., 2019, *Nature*), which establishes proof-of-concept for CAR-T outside oncology but does not derive targets from endometriosis-specific transcriptomic data.

This project asks the question directly from single-cell data: which surface proteins are most specifically expressed in endometriosis lesions relative to healthy endometrial and off-target tissue, at single-cell resolution, and robustly across patients?

We present a differentiable co-expression analysis pipeline built on the Fonseca et al. (2022) endometriosis scRNA-seq atlas that identifies candidate surface targets with substantially better lesion-specificity profiles than the cardiac fibrosis benchmark. For severe, refractory endometriosis — where surgical resection repeatedly fails and quality of life is severely impacted — the risk/reward calculation for an aggressive cellular therapy becomes increasingly compelling, and identifying the right targets is the rate-limiting first step.

The co-expression pair structure of the analysis additionally yields a network-level view of lesion biology, surfacing hub genes around which the endometriosis surface transcriptome is organised.

---

## Key Contributions

- **First scRNA-seq driven surface target identification for endometriosis CAR-T**, using patient-stratified single-cell data rather than bulk RNA or prior literature
- **Better specificity than the cardiac fibrosis benchmark** — top candidates show substantially lower off-target signal than FAP at equivalent enrichment thresholds
- **Co-expression network analysis** as a secondary output — identifying transcriptional hub genes in lesion tissue with potential biological and therapeutic significance

---

## Clinical Context

Endometriosis affects roughly 10% of people with a uterus. In severe (stage III/IV) disease, lesions are difficult to delineate surgically, recurrence rates are high, and hormonal suppression trades efficacy for systemic side effects. For this subset of patients, the risk/reward calculation for an engineered cell therapy shifts — analogous to the tolerance for cytotoxic chemotherapy in oncology where the alternative is worse.

CAR-T targeting of activated fibroblasts in cardiac fibrosis (Aghajanian et al., 2019) establishes the principle that CAR-T can be applied in fibrotic, non-oncological disease with acceptable safety profiles. Our analysis identifies endometriosis-specific targets with a lower estimated off-target footprint than FAP, potentially improving the therapeutic window for this approach.

---

## Method

### Pipeline

```
Raw scRNA-seq (Fonseca et al. 2022, HCA)
        ↓
  Immune cell filtering (CellTypist + Human Endometrium Atlas)
        ↓
  Batch correction & denoising (scVI)
        ↓
  Binarisation (> 1 normalised count per 10k)
        ↓
  Patient-stratified train/test split
        ↓
  Differentiable co-expression energy minimisation
        ↓
  Top surface target extraction + co-expression network analysis
        ↓
  Off-target profiling (Tabula Sapiens multi-organ atlas)
```

### Differentiable Co-expression Scoring

Given a binarised single-cell expression matrix **X** ∈ {0,1}^(N×G) with lesion/healthy labels, we learn soft gene selection weights via a Gumbel-Softmax relaxation:

```
m_j = σ((α_j + g_j) / τ)
```

Per-cell co-activation is computed in log-space:

```
log a_i = Σ_j (1 - X_ij) · log(1 - m_j)
```

The objective maximises lesion activation and minimises healthy activation under a cardinality constraint:

```
L = -E_lesion[log a_i + ε] + E_healthy[a_i] + λ · (Σ_j m_j - k)²
```

The pair search serves two purposes: (1) identifying the single best-specificity target, and (2) mapping co-expression structure across the lesion surfaceome for network analysis.

The loss landscape around the optimum is non-convex but has shallow local minima — the optimiser converges consistently across random seeds, recovering stable high-specificity targets.

---

## Results

Top candidate surface targets on the Fonseca et al. (2022) atlas (patient-stratified, held-out test split):

| Marker | Lesion Enrichment | Notes |
|--------|-------------------|-------|
| GAS1   | High              | Consistently top-ranked across seeds |
| BST2   | High              | Independently emerging as diagnostic biomarker in endometriosis |
| LRP1   | High              | Stable co-expression with GAS1 |
| FAP    | Moderate          | Cardiac fibrosis benchmark — our top candidates exceed this |

Co-expression analysis additionally identifies GAS1, BST2, LRP1, and SGCE as recurring hubs in the lesion surface transcriptome.

---

## Data

No raw data is included in this repository. The analysis uses:

- **Fonseca et al. (2022)** — *A single-cell transcriptomic analysis of endometriosis*, Nature Communications. [Human Cell Atlas](https://explore.data.humancellatlas.org/projects/50db6ba4-3986-4d55-86b7-e1a5a888a17b)
- **Human Endometrium Atlas** — Cell type reference labels. [Cambridge Repository](https://www.repository.cam.ac.uk/handle/1810/373599)
- **Tabula Sapiens** — Multi-organ off-target safety profiling. [tabula-sapiens.sf.czbiohub.org](https://tabula-sapiens.sf.czbiohub.org/)

---

## Installation

```bash
git clone https://github.com/abdurrahman-naeemi/dual-target-cart-endometriosis.git
cd dual-target-cart-endometriosis
pip install -r requirements.txt
```

---

## Usage

```python
from src.model import CoexpressionSelector
from src.data import load_and_preprocess

X, y = load_and_preprocess("path/to/fonseca_denoised.h5ad")

selector = CoexpressionSelector(n_genes=X.shape[1], temperature=0.5, lambda_card=1.0)
selector.fit(X, y, n_epochs=500, lr=0.01)

top_target = selector.get_top_gene()
coexpression_pairs = selector.get_network_pairs(top_k=10)
```

See `notebooks/analysis.ipynb` for the full pipeline.

---

## Repository Structure

```
dual-target-cart-endometriosis/
├── README.md
├── requirements.txt
├── paper/
│   └── geometric_energy_landscapes_cart_endometriosis.pdf
├── src/
│   ├── __init__.py
│   ├── data.py          # scVI denoising, CellTypist filtering, binarisation
│   ├── model.py         # Gumbel-Softmax co-expression selector
│   ├── train.py         # Training loop and seed robustness checks
│   └── evaluate.py      # Specificity scoring, network analysis, off-target profiling
├── notebooks/
│   └── analysis.ipynb
├── figures/
│   └── loss_landscape.png
└── .gitignore
```

---

## Citation

```bibtex
@article{moinudeen2025cart,
  title   = {Geometric Energy Landscapes for Dual-Target CAR-T Marker Discovery in Endometriosis},
  author  = {Moinudeen, Haziq and Naeemi, Abdurrahman},
  year    = {2025},
  note    = {Preprint}
}
```

---

## Authors

**Haziq Moinudeen** — UCL Natural Sciences; Fetal & Placental Physiology Group, EGA Institute for Women's Health  
**Abdurrahman Naeemi** — University of Colorado Boulder (MSc Artificial Intelligence); Lyeons

---

## Acknowledgements

CAR-T in non-oncological fibrotic disease is motivated by Aghajanian et al. (2019). Single-cell data provided by Fonseca et al. (2022) and the Human Cell Atlas. BST2 as an independent endometriosis biomarker: Biomarker Research, 2024.

---

## License

MIT

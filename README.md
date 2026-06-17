# scRNA-seq Driven CAR-T Target Discovery in Endometriosis

Data-driven identification of surface marker candidates for CAR-T therapy in endometriosis, using single-cell transcriptomics, differentiable co-expression analysis, and multi-organ safety profiling against the Tabula Sapiens atlas.

---

## Overview

Most CAR-T target selection outside oncology relies on a small number of candidates proposed from bulk RNA-seq or immunohistochemistry, validated in model organisms, and transferred to humans with limited single-cell resolution. In endometriosis, the closest precedent is FAP targeting from cardiac fibrosis (Aghajanian et al., 2019, *Nature*) — which establishes proof-of-concept for CAR-T in fibrotic non-oncological disease, but does not derive targets from endometriosis-specific transcriptomic data.

This project asks the question directly from single-cell data: which surface proteins are most specifically expressed in endometriosis lesions relative to healthy endometrial and off-target tissue, at single-cell resolution, and robustly across patients?

We present a computational pipeline built on the Fonseca et al. (2022) endometriosis scRNA-seq atlas that identifies candidate surface targets with substantially better lesion-specificity profiles than the cardiac fibrosis benchmark. For severe, refractory (stage III/IV) endometriosis — where hormonal suppression fails, surgical resection recurs, and quality of life is severely impacted — the risk/reward calculation for an aggressive cellular therapy shifts considerably. Identifying the right targets is the rate-limiting first step.

The co-expression structure of the analysis additionally yields a network-level view of the lesion surfaceome, surfacing hub genes around which the endometriosis surface transcriptome is organised — a foundation for graph-theoretic analysis of lesion biology.

---

## Key Contributions

- **First scRNA-seq driven surface target identification for endometriosis CAR-T**, using patient-stratified single-cell data rather than bulk RNA or prior literature
- **Better specificity than the cardiac fibrosis benchmark** — top candidates show substantially lower off-target signal than FAP across 75 tissues in Tabula Sapiens
- **Differentiable co-expression search** over ~4.5 million marker pairs from the 2,647-gene validated surfaceome, recovering biologically stable lesion-enriched combinations
- **Co-expression network structure** as a secondary output — recurring hub genes (GAS1, BST2, LRP1, SGCE) mark stable lesion-associated surface programs and provide a foundation for downstream network analysis

---

## Clinical Context

Endometriosis affects roughly 10% of people with a uterus. In severe (stage III/IV) disease, lesions are difficult to delineate surgically, recurrence rates after excision are high, and hormonal suppression carries significant systemic side effects. For this subset of patients, the risk/reward calculation for an engineered cell therapy is more favourable than it might first appear — analogous to the tolerance for cytotoxic chemotherapy in oncology, where the severity of the disease justifies aggressive treatment.

Aghajanian et al. (2019) establish that CAR-T cells targeting FAP on activated fibroblasts can reduce cardiac fibrosis in mice without detectable off-target toxicity, providing the mechanistic precedent for CAR-T in non-oncological fibrotic disease. Our analysis identifies endometriosis-specific targets with a substantially lower estimated off-target footprint than FAP across human tissues.

---

## Pipeline

```
Bausch-Fluck surfaceome (2,886 proteins)
        ↓
  01_intersect_markers.py
  Intersect with endometriosis atlas gene set
  → 2,647 surface genes present in atlas
  → CSPA mass-spec validated subset
        ↓
  Fonseca et al. (2022) scRNA-seq atlas
  Immune filtering (CellTypist + Human Endometrium Atlas)
  Batch correction + denoising (scVI)
  Patient-stratified binarisation
        ↓
  Differentiable co-expression energy minimisation
  Gumbel-Softmax relaxation over ~4.5M marker pairs
  Gradient descent on log-space AND-gate energy surface
        ↓
  Top surface targets + co-expression network structure
        ↓
  02_downsample_tabula_sapiens.py
  Double-stratified sampling from Tabula Sapiens
  (113,611 cells × 60,606 genes, 75 tissues, RAM-efficient)
        ↓
  03_safety_profile.py
  Patient-stratified co-expression rates across all 75 tissues
  Threshold sweep: 0.2 → 5.0 normalised counts per 10k
  Transfer curves per tissue + heatmap
```

---

## Method

### Surfaceome Intersection

From the Bausch-Fluck 2,886-gene surfaceome, we retain only genes present in the endometriosis atlas and optionally filter to the CSPA mass-spec validated subset (direct experimental evidence of surface localisation). This yields a **2,647-gene candidate space**, giving roughly **4.5 million unordered pairs** to search.

### Differentiable Co-expression Search

For each gene j, a learnable logit α_j encodes selection probability via Gumbel-Softmax:

```
m_j = σ((α_j + g_j) / τ),   m_j ∈ (0, 1)
```

Per-cell co-activation is computed in log-space for numerical stability:

```
log a_i = Σ_j (1 - X_ij) · log(1 - m_j)
```

The objective maximises lesion activation and suppresses healthy activation under a cardinality constraint:

```
L = -E_lesion[log a_i + ε] + E_healthy[a_i] + λ · (Σ_j m_j - k)²
```

This frames marker discovery as geometry on a data-driven energy surface rather than brute-force enumeration.

**Optimisation landscape.** The loss surface is non-convex at the scale of parameter space we are working with, but has shallow local minima — the optimiser converges consistently across random seeds, recovering biologically stable marker combinations.

![Loss Landscape](figures/loss_landscape.png)

### Safety Profiling

Candidate markers are profiled against the Tabula Sapiens multi-organ human cell atlas (113,611 cells, 75 tissues). Co-expression rates are computed with patient stratification — rates are computed per donor within each tissue, then averaged across donors — to prevent single-donor outliers from dominating the signal.

Rather than reporting at a single threshold, we generate **transfer curves**: co-expression rate as a function of expression threshold across the full range 0.2 → 5.0 normalised counts per 10k. This gives a calibration-curve view of the safety profile that does not depend on an arbitrary binarisation choice.

---

## Results

### Top Candidate Pairs (Exact Enumeration)

| Rank | Marker Pair     | Specificity Score |
|------|-----------------|-------------------|
| 1    | GAS1 + LRP1     | 0.579             |
| 2    | BST2 + GAS1     | 0.534             |
| 3    | ATRAID + GAS1   | 0.514             |
| 4    | CD164 + GAS1    | 0.500             |
| 5    | DDR2 + GAS1     | 0.481             |
| 6    | BST2 + LRP1     | 0.477             |
| 7    | IL6ST + SGCE    | 0.470             |
| 8    | PRNP + SGCE     | 0.470             |

GAS1 appears in 5 of the top 8 pairs. BST2, LRP1, and SGCE recur across the top 10. This consistency across pairs — rather than a single dominant marker — supports the view that these genes mark a stable lesion-associated surface program in this dataset.

The differentiable approach recovers GAS1 + BST2 with approximately **10× lesion enrichment** in a held-out test split. GAS1 + BST2 is the second-best exact pair, suggesting the energy landscape has shallow local minima that still guide the search toward highly enriched targets.

BST2 is independently emerging as a diagnostic biomarker in endometriosis (Biomarker Research, 2024), providing external validation that our approach is recovering biologically meaningful signal.

### Multi-Organ Safety Profile

Patient-stratified co-expression rates for GAS1 + BST2 at primary threshold >1.0 normalised counts per 10k, across 75 human tissues:

```
Tissue                                    Co-expr%   Flag
--------------------------------------------------------------
right ovary                                ~4-15%    WARN (threshold-dependent)
left ovary                                 ~4-13%    WARN (threshold-dependent)
adipose tissue                             ~3-11%    WARN (threshold-dependent)
heart right ventricle                      <7%       SAFE at threshold >1.0
endometrium                                0.19%     SAFE
uterus                                     0.00%     SAFE
blood, bone marrow, spleen                 <0.2%     SAFE
kidney, liver, lung                        <0.4%     SAFE
```

The ovarian signal is highest, which is biologically expected given anatomical and developmental proximity to endometriosis. At the primary threshold (>1.0), critical organs (heart, kidney, liver, lung) are well within safe range. The full threshold sweep is shown in the transfer curves below.

![Safety Profile — Patient-Stratified Transfer Curves](figures/transfer_curves.png)

---

## Repository Structure

```
endometriosis-cart-target-discovery/
│
├── README.md
├── requirements.txt
├── .gitignore
│
├── paper/
│   └── geometric_energy_landscapes_cart_endometriosis.pdf
│
├── src/
│   ├── 01_intersect_markers.py          # surfaceome × atlas intersection → marker shortlist
│   ├── 02_downsample_tabula_sapiens.py  # RAM-efficient double-stratified atlas downsampling
│   └── 03_safety_profile.py             # patient-stratified safety profiling + transfer curves
│
├── data/
│   ├── surfaceome_in_atlas.csv          # 2,647 surface genes present in endometriosis atlas
│   ├── surfaceome_in_atlas_CSPA_validated.csv   # mass-spec validated subset
│   └── marker_shortlist.txt             # plain gene list
│
└── figures/
    ├── transfer_curves.png              # patient-stratified safety profile across 75 tissues
    └── loss_landscape.png              # optimisation energy landscape around optimal pair
```

---

## Installation

```bash
git clone https://github.com/abdurrahman-naeemi/endometriosis-cart-target-discovery.git
cd endometriosis-cart-target-discovery
pip install -r requirements.txt
```

---

## Data

No raw data is included in this repository. The analysis uses:

- **Fonseca et al. (2022)** — *A single-cell transcriptomic analysis of endometriosis*, Nature Communications. [Human Cell Atlas](https://explore.data.humancellatlas.org/projects/50db6ba4-3986-4d55-86b7-e1a5a888a17b)
- **Human Endometrium Atlas** — Cell type reference labels. [Cambridge Repository](https://www.repository.cam.ac.uk/handle/1810/373599)
- **Bausch-Fluck surfaceome** — 2,886-gene validated human surfaceome. [surfaceome.org](http://wlab.ethz.ch/surfaceome/)
- **Tabula Sapiens** — 113,611 cells across 75 human tissues. [tabula-sapiens.sf.czbiohub.org](https://tabula-sapiens.sf.czbiohub.org/)

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

CAR-T in non-oncological fibrotic disease: Aghajanian et al. (2019), *Nature*. Single-cell endometriosis data: Fonseca et al. (2022), *Nature Communications*. BST2 as independent endometriosis biomarker: *Biomarker Research*, 2024. Gumbel-Softmax relaxation: Jang, Gu & Poole, ICLR 2017. scVI denoising: Lopez et al., *Nature Methods* 2018.

---

## License

MIT

# GCT-AVP framework  
**A Generative–Constraint–Targeted Framework for Orally Active Antiviral Peptide Discovery**
---

## Overview

GCT-AVP is a **hierarchical multi-objective computational framework** for de novo discovery of orally active antiviral peptides (AVPs). It integrates **generative sequence exploration**, **biophysical constraint filtering**, and **target-specific structural prioritization** to solve the core translation gap in peptide drug design: balancing antiviral potency, host safety, and physiological stability.

---

## Framework architecture

Fully aligned with the manuscript, the framework sequentially narrows large peptide sequence spaces into high-confidence candidates:

### 1. Generative Exploration
- MSA-informed diffusion model (Evodiff-MSA_OA_DM_MAXSUB)
- Large-scale de novo sequence generation (15–25 amino acids)
- Evolutionary conservation-guided sequence sampling

### 2. Constraint-Driven Selection
- Physicochemical filtering (net charge, hydrophobicity, stability, amphipathicity)
- Antiviral activity prediction (Stack-AVP)
- Cytotoxicity risk elimination (ToxinPred3.0)
- Multi-constraint sequence space compression

### 3. Target-Specific Evaluation
- PEDV Spike protein-targeted activity scoring (LSTM regression)
- De novo peptide structure prediction (PEP-FOLD4)
- Molecular docking & binding stability assessment (HDOCK)
- Multi-objective candidate prioritization

---

## Repository structure
```
GCT-AVP-framework/
├── data/ 
│ ├── example_sequences.csv # Sample peptide sequences
│ └── processed_features.csv # Precomputed physicochemical features
│
├── pipeline/ 
│ ├── step1_preprocess.py 
│ ├── step2_physicochemical_filter.py 
│ ├── step3_feature_extraction.py 
│ └── step4_ranking.py 
│
├── generative/ 
│ └── sequence_generation.py 
│
├── selection/ 
│ ├── constraint_filtering.py 
│ └── cytotoxicity_evaluation.py 
│
├── evaluation/ 
│ ├── target_scoring.py 
│ └── multi_objective_ranking.py 
│
├── results/ 
│ ├── filtered_sequences.csv 
│ └── prioritized_candidates.csv 
│
├── docs/ 
│ └── reproducibility.md 
│
├── environment.yml 
├── requirements.txt 
├── LICENSE 
└── README.md 
```

---

## Getting Started

### Installation Guide

Welcome to the GCT-AVP framework. This guide describes how to set up the environment and run the reproducible components of the pipeline.

### Prerequisites
- Python 3.11
- Conda (recommended for environment management)
- CPU/GPU (GPU recommended for generative modeling)

Required libraries:`numpy`, `pandas`, `tqdm`, `scikit-learn`, `xgboost`.

### Setting Up the Environment

1. **Clone the repository**

```bash
git clone https://github.com/HanranZhao-spec/GCT-AVP-framework.git
cd GCT-AVP-framework
 ```

2. **Create a Conda Environment** 
 ```bash
conda create -n gct-avp python=3.11
conda activate gct-avp
   ```

3. **Install dependencies**  
 ```bash
pip install -r requirements.txt
   ```

4. **Install PyTorch and Related Packages**  
EvoDiff requires PyTorch along with additional libraries. The following example demonstrates the installation of a CPU-compatible version of PyTorch. For optimal performance, adjust the pytorch version based on your system’s specifications. Install the required packages using the following commands:
   ```bash
   conda install pytorch torchvision torchaudio cpuonly -c pytorch
   ```

---

## Pipeline execution

Run the full GCT-AVP filtering pipeline sequentially:

### Input

A CSV file containing peptide sequences:
CSV file with peptide sequences (data/example_sequences.csv)

---

### Step 1: Sequence preprocessing

```bash
python pipeline/step1_preprocess.py
 ```

### Step 2: Physicochemical filtering
```bash
python pipeline/step2_physicochemical_filter.py
 ```

### Step 3: Feature extraction
```bash
python pipeline/step3_feature_extraction.py
 ```

### Step 4: Candidate ranking
```bash
python pipeline/step4_ranking.py
 ```

### Output

Top-ranked candidates (results/prioritized_candidates.csv)

---

## Framework logic

The GCT-AVP framework implements a **constraint-based filtering strategy** rather than a single predictive model.

Key principles:
- Sequential reduction of sequence space  
- Multi-objective optimization (activity, safety, stability)  
- Integration of data-driven and knowledge-based constraints  

This design enables efficient prioritization of experimentally testable candidates from large generated libraries.

## External tools

These components were accessed via publicly available servers or standalone implementations and are described in detail in the manuscript.
- Antiviral prediction models (Stack-AVP)  
- Cytotoxicity prediction tools (ToxinPred3.0)  
- Structure prediction (PEP-FOLD4)  
- Molecular docking (HDOCK server)  

This repository contains the self-contained core pipeline independent of external server dependencies.

---

## Reproducibility

All scripts and parameters are fully reproducible. Detailed step-by-step instructions for reproducing the manuscript’s computational results are provided in:
docs/reproducibility.md

## Design philosophy

GCT-AVP emphasizes:

- integration of heterogeneous computational tools  
- constraint-driven candidate prioritization  
- experimental validation as the ultimate selection criterion  

This framework bridges large-scale sequence generation and practical antiviral peptide discovery.

---

## Contributing

Contributions are welcome! Please open an issue or submit a pull request for any improvements or bug fixes.

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

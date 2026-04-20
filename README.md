# GCT-AVP framework  
A Generative–Constraint–Targeted framework for antiviral peptide discovery

---

## Overview

GCT-AVP is a hierarchical computational framework for antiviral peptide (AVP) discovery, integrating generative exploration, multi-constraint virtual screening, and target-specific prioritization.

The framework is designed to efficiently navigate large peptide sequence spaces while enforcing biologically relevant constraints, enabling the identification of candidate peptides with favorable antiviral activity and safety profiles.

This repository provides a **reproducible implementation of the core computational pipeline** described in our study.

---

## Correspondence to manuscript

This repository reproduces key computational components reported in the manuscript, including:

- Sequence preprocessing and filtering  
- Physicochemical characterization  
- Multi-parameter candidate ranking  
- Candidate diversity analysis  

Experimental validation and external prediction tools are described in the manuscript.

---

## Framework architecture

The GCT-AVP framework consists of three modules:

### 1. Generative exploration
- Evodiff-MSA_OA_DM_MAXSUB
- Large-scale sequence space exploration

### 2. Constraint-driven selection
- Physicochemical filtering (charge, hydrophobicity, stability)
- Functional prediction (external tools)
- Cytotoxicity evaluation (external tools)

### 3. Target-specific evaluation
- Target-aware scoring
- Structure-based assessment (external tools)
- Multi-objective prioritization

---

## Repository structure
```
GCT-AVP-framework/
├── data/
│ ├── example_sequences.csv
│ └── processed_features.csv
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
├── environment.yml
├── requirements.txt
├── LICENSE
└── README.md
```

## Getting Started

### Installation Guide

Welcome to the GCT-AVP framework. This guide describes how to set up the environment and run the reproducible components of the pipeline.

### Prerequisites

- Python 3.11  
- Conda (recommended)  

Required libraries:`numpy`, `pandas`, `tqdm`, `scikit-learn`, `xgboost`.

### Setting Up the Environment

1. **Clone the repository**

```bash
git clone https://github.com/HanranZhao-spec/GCT-AVP-framework.git
cd GCT-AVP-framework
 ```

2. **Create a Conda Environment** 
 ```bash
conda create --name gct-avp python=3.11
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

The GCT-AVP framework operates as a sequential filtering and prioritization pipeline.

### Input

A CSV file containing peptide sequences:

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

These tools were accessed via their official web servers or standalone implementations as described in the manuscript.

This repository focuses on the **reproducible core pipeline** independent of these external dependencies.

---

## Reproducibility

All scripts provided in this repository are sufficient to reproduce the computational filtering and ranking workflow described in the study.

Detailed step-by-step instructions are available in:
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

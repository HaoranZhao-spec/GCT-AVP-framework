# GCT‑AVP framework
**A Generative–Constraint–Targeted Framework for Broad‑Spectrum Antiviral Peptide Discovery**

---

## Overview
GCT‑AVP is a hierarchical, multi‑objective computational framework for de‑novo discovery of broad‑spectrum orally‑active antiviral peptides (AVPs).
This framework addresses key bottleneck in peptide therapeutic development: simultaneous optimization for broad‑spectrum antiviral activity, host safety and physiological stability.

By integrating generative modeling, constraint‑driven filtering and broad‑spectrum enveloped‑virus targeted predictive scoring, GCT‑AVP navigates huge peptide sequence space to deliver experimentally tractable antiviral peptide candidates.

>  Important note:
>  **Locally runnable**: MSA‑diffusion generation, physicochemical calculation, composite potency scoring
>  **External web‑server only**: Antiviral prediction (Stack‑AVP), cytotoxicity prediction (ToxinPred3.0). These two steps cannot run offline locally, please follow protocol files.

## Framework architecture
Sequential constraint‑based multi‑step pipeline:
### 1. Generative Exploration
- MSA‑informed diffusion model (`Evodiff‑MSA_OA_DM_MAXSUB`)
- Large‑scale de‑novo peptide generation (15–25 amino acids length)
- Evolutionary‑conservation guided sequence sampling

### 2. Constraint‑Driven Selection
1. Physicochemical filtering: net charge, hydrophobicity, aliphatic index, instability index, isoelectric point etc.
2. Antiviral activity prediction via Stack‑AVP web server
3. Cytotoxicity risk elimination via ToxinPred3.0 web server

### 3. Target‑Specific Evaluation
Multi‑strain enveloped virus potency scoring (ESM‑2‑Ridge derived composite score), output ranked candidate list.

## Repository structure```
GCT-AVP-framework/
│
├── pipeline/ 
│ ├── step1_Generative_exploration 
│ ├── step2_Physicochemical_constraint
│ ├── step3_Functional_constraint 
│ ├── step4_Cytotoxicity_constraint 
│ └── step5_Target_specific_potency_score
│
├── Generative_Exploration/ 
│ └── msa_diffusion_generation.py 
│
├── Constraint_driven_Selection/ 
│ ├── Physicochemical filtering.py 
│ ├── Functional constraint
│        ├── input_sequences.csv
│        ├── output_predictions.csv
│        └──  protocol.md
│ ├── Cytotoxicity constraint
│        ├── input_sequences.csv
│        ├── output_toxicity.csv
│        └──  protocol.md
│
├── Target specific potency score/ 
│        ├──Specific Potency scoring.py 
│        └── supplement figs.py
│ 
├── results/ 
│ └──  final candidate.csv 
│
├── docs/ 
│ └── reproducibility.md 
│
├── environment.yml 
├── requirements.txt 
├── LICENSE 
├── CITATION.cff 
├── .gitignore
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

---

### Step 1: Generative exploration
Input: directory with .a3m MSA multiple sequence alignment files

```bash
python generative/msa_diffusion_generation.py \
  --directory_path ./data/example_msa \
  --output_csv_file ./generated_peptides.csv \
  --to_device cuda \
  --total_sequences 200
 ```

### Step 2:  Physicochemical property calculation
```bash
python physicochemical_filtering.py \
  --input_csv_file ./generated_peptides.csv \
  --output_csv_file ./physicochem_output.csv
 ```

### Step 3: Antiviral activity classification (WEB‑SERVER ONLY)
Export sequences from step2 output CSV as input_sequences.csv into selection/Functional_constraint/
Follow instructions: selection/Functional_constraint/protocol.md
Upload file to Stack‑AVP web‑server. Download prediction output. Keep sequences predicted as AVP.

### Step 4: Cytotoxicity prediction (WEB‑SERVER ONLY)
Filter AVP‑positive sequences from step3, export to selection/Cytotoxicity_constraint/input_sequences.csv
Follow: selection/Cytotoxicity_constraint/protocol.md
Retain only sequences predicted as non‑toxic.

### Step 5: Target‑oriented composite potency scoring 
```bash
python Target_specific_potency_score/specific_potency_scoring.py \
  --input ./scoring_results_fixed.csv \
  --output ./results/prioritized_candidates.csv
 ```

Final ranked candidates saved in results/prioritized_candidates.csv, sorted descending by composite score.

---

## Framework logic

GCT‑AVP uses constraint‑based sequential filtering strategy, not single end‑to‑end predictive model:

Key principles:
- Sequential reduction of sequence space  
- Multi-objective optimization (activity, safety, stability)  
- Combine knowledge‑based constraints and data‑driven machine‑learning predictions  

This design enables efficient prioritization of experimentally testable candidates from large generated libraries.

## External tools

These components were accessed via publicly available servers or standalone implementations and are described in detail in the manuscript.
- Antiviral prediction models (Stack-AVP)  
- Cytotoxicity prediction tools (ToxinPred3.0)  


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

\# Reproducibility for GCT‑AVP framework

Full computational reproduction workflow for manuscript results.



\## Hardware requirement

GPU: ≥10GB VRAM (recommended). CPU mode supported but slow for sequence generation.

Python 3.11.



\## Full workflow steps

1\. Run generative MSA‑diffusion peptide generation

&#x20;   `python generative/msa\_diffusion\_generation.py ...`

&#x20;   Output: generated peptide csv.



2\. Compute physicochemical descriptors

&#x20;   `python physicochemical\_filtering.py ...`



3\. Stack‑AVP antiviral prediction (web‑server only)

Follow `selection/Functional\_constraint/protocol.md`.

Filter sequences predicted as AVP.



4\. ToxinPred3.0 cytotoxicity filtering (web‑server only)

Follow `selection/Cytotoxicity\_constraint/protocol.md`.

Retain only non‑toxic entries.



5\. Target‑specific composite potency scoring

Run specific\_potency\_scoring.py.

Final ranked candidates output to results/prioritized\_candidates.csv.



\## Note

Step3 and Step4 depend external public web servers, cannot run offline locally.

Minor score fluctuation may appear due server side model update.

All parameters used in manuscript are default as described in protocol.md.


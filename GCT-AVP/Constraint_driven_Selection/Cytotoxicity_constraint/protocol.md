ToxinPred3.0 Cytotoxicity Prediction Protocol

Tool Information

Tool: ToxinPred3.0

Access: Public web server

Purpose: Prediction of peptide toxicity

Input Data



Input sequences are provided in:



input\_sequences.csv



Format:



One peptide per row

Column name: sequence

Procedure

Access the ToxinPred3.0 web server

Upload or paste peptide sequences

Select toxicity prediction mode

Use default parameters unless otherwise specified

Run prediction

Export results as CSV

Output



Results are provided in:



output\_toxicity.csv



Typical fields:



sequence

toxicity prediction (toxic / non-toxic)

confidence score (if available)

Filtering Criteria

Sequences predicted as "toxic" were excluded

Only "non-toxic" candidates were retained for downstream analysis

Reproducibility Notes

Default model and parameters were used

Output files provided here correspond exactly to those used in the manuscript

Limitations



As a web-based tool, model updates may affect reproducibility over time.

All results here reflect the state of the server at the time of analysis.


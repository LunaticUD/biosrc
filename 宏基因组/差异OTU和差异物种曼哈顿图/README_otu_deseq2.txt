OTU differential-abundance and Manhattan-plot example

Files
-----
- otu_deseq2_manhattan.py: complete Python conversion
- OTU_example.txt: example OTU count table with taxonomy in the final column
- metadata_example.txt: example metadata

Install recommended backend
---------------------------
pip install pydeseq2 pandas numpy scipy matplotlib statsmodels

Run
---
python otu_deseq2_manhattan.py

Using your own files
--------------------
Set:
WORK_DIR = Path(r"C:/Users/14266/Desktop/OTU")
OTU_FILE = WORK_DIR / "OTU.txt"
METADATA_FILE = WORK_DIR / "metadata.txt"

The metadata file must contain:
sample    Group

Comparison direction
--------------------
The script compares:
CON Day 28 / CON Day 0

Therefore:
- positive log2FoldChange: higher on CON Day 28
- negative log2FoldChange: higher on CON Day 0

Backend
-------
The script uses PyDESeq2 when installed. If it is unavailable, the default
'auto' setting uses a negative-binomial GLM approximation for testing the full
workflow. Install PyDESeq2 for the intended DESeq2-style analysis.

Font
----
All plot text is configured as Times New Roman.

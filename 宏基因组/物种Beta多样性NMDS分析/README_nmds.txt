Bray-Curtis NMDS analysis

Files:
- nmds_bray_analysis.py: complete Python conversion
- tax.txt: example taxa-by-sample abundance table
- group.txt: example sample grouping table

Install dependencies:
pip install pandas numpy scipy matplotlib scikit-learn

Run:
1. Put nmds_bray_analysis.py, tax.txt and group.txt in the same folder.
2. Open a terminal in that folder.
3. Run:
   python nmds_bray_analysis.py

Outputs:
- bray-crutis.txt
- pcoa_coordinates_nmds_workflow.csv
- nmds_coordinates.csv
- nmds_summary.csv
- NMDS_stressplot.png
- NMDS1.png
- NMDS2.png

Font:
All plot text is configured as Times New Roman.
On Windows, Times New Roman is normally available by default.

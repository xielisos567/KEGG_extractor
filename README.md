#Note: The annotated files by KEGG tools such as KofamKOALA are needed for gene extraction. The files should be processed like the format in the example document.

# Usage of KEGG_extractor:
python3  KEGG_extractor.py  -i  example/protein  -f  example/Wood-Ljungdahl_20.txt  -s example/species  -o  example/results/result_protein_WL

python3  KEGG_extractor.py  -i example/cds  -f  example/Wood-Ljungdahl_20.txt  -s example/species  -o  example/results/result_cds_WL

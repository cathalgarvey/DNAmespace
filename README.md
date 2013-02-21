# DNAmespace
A Python module for presenting bacterial genomes (from NCBI/Genbank) as namespaces in Python.

Usage (in iPython):
In [0]: import dnamespace
In [1]: ecoli_w3110 = dnamespace.new("test_genomes/E.coli_K12_W3110.gbk")
In [2]: ecoliw3110.lac<tab>
ecoliw3110.lacA  ecoliw3110.lacI  ecoliw3110.lacY  ecoliw3110.lacZ


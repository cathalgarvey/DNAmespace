# DNAmespace
Copyright: Cathal Garvey
License: GNU Affero General Public License v3
(License text included as dnamespace.license in Python prompt)

## What
A Python module for presenting bacterial genomes (from NCBI/Genbank) as
namespaces in Python, allowing tab-completion of genes and other named features.

## Why
There are many tools for studying wild DNA, and there are many for manipulating
DNA in a pleasant graphical environment, but there are few in-betweens for
people who like to work with logical code structures. This is intended to be
a building block for people who would rather a "bioprogramming language",
allowing them to access the vast databases of bacterial genomes available
online and to use them in Python as namespace objects for data lookup.

For example, if one wished to "refactor" a bacterial genome, one could import
the wild-type genome of the target species, call all of the desired genes,
regulatory regions and other DNA features they wish to refactor, and reassemble
them into a target genome, in a prompt environment. Perhaps it would look like
this:
> import genomebuilder # Not yet implemeted
> import dnamespace
> NewGenome = genomebuilder.new(codon_table="11") # Not yet implemented
> with dnamespace.new("oldgenome") as Source: # Not yet implemented
>     for gene in Source._gene_list: # Not yet implemented.
>         NewGenome.add_gene(gene) # Not yet implemented
>     for ncRNA in Source._ncrna_list: # Not yet implemented.
>         NewGenome.add_ncRNA(ncRNA) # Not yet implemented
>     # And all the other stuff one might copy over..
>     NewGenome.ori = Source.origin # Not yet implemented
> NewGenome.add_gene(genomebuilder.viral_bootstrapper(NewGenome)) # Not yet implemented
> NewGenome.add_signature("This genome was built by Cathal.") # Not yet implemented
> with open("newgenome",mode="w") as Output:
>     Output.write(str(NewGenome)) # Not yet implemented

## Usage
(Below sample in/out is in iPython, tab completion may not work in your
preferred prompt)

In [0]: import dnamespace

In [1]: ecoli_w3110 = dnamespace.new("test_genomes/E.coli_K12_W3110.gbk")

DNAmespace objects present genome features (presently only genes) as
object attributes, enabling tab-completion in iPython and other
smart shells:

*In [2]: ecoliw3110.lac<tab>*

_ecoliw3110.lacA  ecoliw3110.lacI  ecoliw3110.lacY  ecoliw3110.lacZ_

*In [3]: ecoliw3110.lacZ.<tab>*

_Ecoli.lacZ.amino        Ecoli.lacZ.meta         Ecoli.lacZ.sequence_
_Ecoli.lacZ.aminos       Ecoli.lacZ.orfs         Ecoli.lacZ.transcript_
_Ecoli.lacZ.features     Ecoli.lacZ.rna          Ecoli.lacZ.transcripts_


To ensure compatibility with operons as well as singular genes,
gene objects retain lists for transcripts (CDSs and soon RNA),
orfs (protein-coding subset of transcripts), rnas (non-protein-coding
rnas), and translated CDSs. For each list (designated by the plural),
there is a singular property which returns the first item in the list,
or None if the list is empty.

*In [4]: ecoliw3110.lacZ.transcripts*

*Out[4]: ['ATGACCATGATTACGGATTCACTG...TGTCAAAAATAA']*

*In [5]: ecoliw3110.lacZ.transcript*

*Out[5]: 'ATGACCATGATTACGGATTCACTG...TGTCAAAAATAA'*


## Code Layout
You might want to re-use bits of DNAmespace separately, and I'll try to
keep code dependencies sane rather than webbed so that useful sub-units
can be re-used. Perhaps you want to namespace-ify something unrelated?
Easy, just copy virtualns and gnulicenses into your project and subclass
virtualns.nsdict.

* nucutils and virtualns are both standalone, but require gnulicenses.
* Import of genbank files is managed (hideously) by parsegb, which requires nucutils and gnulicenses.
* genomespace requires nucutils and parsegb, in addition to gnulicenses.

## Todo
DNAmespace isn't even remotely finished.
* At present, DNAmespace only parses genes, and of those genes it only imports "CDS" features, not RNAs.
* The parsegb module is spaghetti-code and needs refactoring.
* The core genomespace.genomespace object should be a virtualns.nsdict subclass.

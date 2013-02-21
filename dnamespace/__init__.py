'''DNAmespace, a namespace-like interface to bacterial genomes.
by Cathal Garvey, licensed under the GNU Affero GPLv3 license.
License text can be accessed via dnamespace.license
Email: cathalgarvey@cathalgarvey.me
Twitter: @onetruecathal
Code: https://gitorious.org/~cathalgarvey
Blog: http://www.indiebiotech.com

DNAmespace is part of a general push towards creating a pythonic toolset
for designing DNA constructs at a scale between vectors and whole genomes.
The intention of DNAmespace is to provide an object that, given a genbank
file such as a bacterial chromosome/genome downloaded from NCBI, will
present an interface akin to an imported library, so that users can
access features and genes by name, for example:
>>> import dnamespace
>>> ecoli = dnamespace.new("E.coli_K12_W3110.gbk") # Importing from file
>>> laczCDS = ecoli.lacz.cds # CDS access
>>> laczaminos = ecoli.lacz.translation # CDS translation
>>> my_foo_sequence = ecoli[40000:78000] # Direct genome addressing
>>> my_foo_reverse_complement = ecoli[40000:78000].reverse_complement()
>>> rhodopsin_cds = ecoli.rhodopsin.cds


To account for widely varying differences in the way genbank files are
marked up for genes, whether single-protein or polycistronic, the
DNAmespace "api" for genes assumes that all genes are polycistronic,
and has a list property called "cistrons". The "cds" property in the
case of any given gene, polycistronic or no, contains the *first* CDS
feature in the *feature table*, not necessarily the first in-sequence.
'''
from dnamespace import genomespace
from dnamespace.gnulicenses import Affero as license

def new(filen):
    return genomespace.genomespace(filen)

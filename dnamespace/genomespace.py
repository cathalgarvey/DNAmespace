#!/usr/bin/env python3
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
from dnamespace import parsegb
from dnamespace import virtualns
from dnamespace.gnulicenses import Affero as license

# ecoli.<tab>
# ecoli.geneN - All gene names represented as attributes.
# ecoli.geneN_ - Gene names that collide with keywords are suffixed with "_"
# ecoli.fooB.<tab>
# transcript    - the first transcript, if any.
# transcripts   - a subclassed list of all CDSs.
# amino         - the first transcript's translation.
# aminos        - a subclassed list of translations.
# meta          - Metadata from genbank feature entry.
# __doc__       - Set to one of the key feature table meta descriptors, like "note"
import keyword

class geneNS(virtualns.nsdict):
    'Present a namespace or dict interface for a gene.'
    # Only "db_xref; locus_tag; gene" were common to all genes in E.coli DH10B.
    # These keys are not necessarily common to all sub-gene features, and are
    # not common to *all* features (i.e. "gene" will not be found in "origin")
    def __init__(self, *args, **nargs):
        # autofix means any keys that conflict with methods or keywords
        # will be suffixed with underscores until they don't. This means
        # that, for example, the e.coli gene "def" will become "def_".
        virtualns.nsdict.__init__(self, autofix=True)

        # Common containers / Gene API:
        # sequence will contain the whole-gene sequence, probably from
        # the "gene" feature.
        self['sequence'] = ''
        # transcripts are rna copies of DNA, and are either from CDS
        # or RNA features.
        self['transcripts'] = []
        # orfs and rnas are subsets of transcripts: Entries in either list
        # will always be duplicated in transcripts.
        self['orfs'] = []
        self['rna'] = []
        self['features'] = []
        self['aminos'] = []
        self['meta'] = {}

    @property
    def transcript(self):
        if self['transcripts']:
            return self['transcripts'][0]

    @property
    def amino(self):
        if self['aminos']:
            return self['aminos'][0]

    def _import_gbfeature(self, gbfeature):
        'Enters a GBFeature object into namespace structure.'
        # Feature types in testgenomes:
        # Feature types encountered across test genomes:
        # source; gene; CDS; mobile_element; ncRNA; misc_feature; STS;
        # misc_RNA; rRNA; tRNA; repeat_region; tmRNA; mat_peptide;
        # rep_origin; exon; polyA_signal; polyA_site

        # Determine Feature Type
        if gbfeature.type == "gene":
            # Initiate and store gene feature meta.
            self._handle_gene(gbfeature)
        elif gbfeature.type == "CDS":
            # Add to transcripts list.
            self._handle_CDS(gbfeature)
        elif gbfeature.type == "misc_feature":
            # Add to misc_features list.
            self._handle_misc_feature(gbfeature)
        else:
            # Deal with as "other".
            self._handle_other_feature(gbfeature)

    def _handle_CDS(self, gbfeature):
        'Imports CDS-specific data.'
        # Should possibly try to determine order relative to existing
        # entries, *without* calling the sequence property (to save CPU/RAM).
        # Might need to refactor parsegb to allow intermediate range calcs:
        # i.e. features are parsed to get range in a native format at init,
        # but sequence isn't read until called. Range can then easily be
        # read/compared by client software like dnamespace.
        self['transcripts'].append(gbfeature.sequence)
        self['orfs'].append(gbfeature.sequence)
        # Need to handle translations and put in self['aminos'].
        try:
            self['aminos'].append(gbfeature.meta['translation'])
        except KeyError:
            NotImplementedError("No translation found for CDS {0}, and translation from nucleotides is not yet implemented.".format(gbfeature.meta['gene']))

    def _handle_gene(self, gbfeature):
        'Imports gene-specific data.'
        # Need to import gbfeature.meta to self['meta'].
        pass

    def _handle_misc_feature(self, gbfeature):
        'Imports misc_feature specific data.'
        self['features'].append(gbfeature)

    def _handle_other_feature(self, gbfeature):
        'Imports data from a non-gene/CDS/misc_feature feature.'
        # This is a temporary hack and *will* change. Features is for
        # misc_features, and noncoding RNA belongs in transcripts.
        self['features'].append(gbfeature)

    def _finalise(self):
        '''Called after all related Features are grouped in this namespace.
        Seeks out critical data that applies to the gene and restructures
        it accordingly; the first transcripts entry becomes transcript,
        first aminos entry becomes amino, and most-relevant-docstring
        becomes self.__doc__.'''
        # Get first transcript and amino (if any), where functinoal RNA
        # counts as a transcript:
        pass
        # Find most-significant-meta for __doc__.
        # Meta types found in test_genomes set:
        # sub_strain; strain; mol_type; organism; db_xref; locus_tag; gene;
        # codon_start; product; transl_table; note; translation; protein_id;
        # mobile_element_type; gene_synonym; ncRNA_class; pseudo; allele;
        # standard_name; function; rpt_type; rpt_unit_range; EC_number;
        # transl_except; ribosomal_slippage; GO_process; GO_component;
        # experiment; GO_function; old_locus_tag; codon_recognized;
        # pseudogene; inference; NonD; TIGR00976; chromosome; map; number;
        # collection_date; country; isolation_source; rpt_family
        # ---
        # Probably useful metadata:
        # note: Created by GENtle
        # function: self-explanatory
        # product: Often describes the translation product.

class genomespace:
    '''Provides a namespace-like object interface to a genbank file.'''
    def __init__(self, gb_file, keepfile=False):
        'Should be created with a path or filename for a valid genbank file.'
        self._gbfile = parsegb.GenbankFile(file_name=gb_file)
        # GenbankFiles have a .features list containing GBFeature objects
        # The GBFeature meta dict will usually contain a "gene" key:
        # Actual gene entries have this, and sub-parts of the gene will
        # usually do so, too. So, DNAmespace first parses features to
        # subordinate anything with a "gene" entry into a dict keyed by
        # gene name.
        self._subordinate_genes()
        # Make all gene dict entry names into properties:
        self._make_gene_properties()
        # If sequences have already been looked-up and internalised,
        # then we can save some RAM by deleting the original gbfile
        # object.
        if not keepfile:
            del(self._gbfile)

    def _subordinate_genes(self):
        'Parse genbank file features and organise by gene name.'
        self._genes = {}
        for feature in self._gbfile.features:
            try:
                gene_name = feature.meta['gene']
                # If no gene name is found, none of the below happen..
                if gene_name not in self._genes.keys():
                    # Create new geneNS object in self._genes.
                    # self._genes keys are not directly exposed, they
                    # are first parsed to suffix reserved keywords with "_",
                    # so a gene might be called "def" in self._genes, but
                    # be self.def_ when called as a property.
                    self._genes[gene_name] = geneNS()
                # Now that we know there's a geneNS object with this name
                # we'll pass the feature to the geneNS object's handler.
                self._genes[gene_name]._import_gbfeature(feature)

            except KeyError:
                # For now, we don't care about features with incomplete
                # metadata and no "gene" meta; important or useful genes
                # are likely to have correct meta.
                # Would be nice to write a "last resort" method to find
                # requested features that weren't successfully extracted
                # if possible? Perhaps appending this method to a custom
                # __getattr__ method and making it contingent on an
                # instance argument, like "desperate=True"?
                pass

    def _make_gene_properties(self):
        for gene_name in self._genes.keys():
            if gene_name in keyword.kwlist:
                # This is ugly but necessary, sadly..
                new_gene_name = gene_name + "_"
            else:
                new_gene_name = gene_name
            # Make this less simplistic; features must have a namespace-
            # like interface, too.
            self.__dict__[new_gene_name] = self._genes[gene_name]

    def _make_feature_namespace(self, gene_name):
        'Presents a namespace interface to a self.genes dict entry.'
        # This returning of dict objects is a dumb space-filler for now.
        return self._genes[gene_name]

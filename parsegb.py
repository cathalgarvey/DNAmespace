#!/usr/bin/env python3
'''parsegb - A collection of objects for handling genbank files.
by Cathal Garvey;
Email: cathalgarvey@cathalgarvey.me
Twitter: @onetruecathal
Code: https://gitorious.org/~cathalgarvey
Blog: http://www.indiebiotech.com

Parsegb is designed around one core class, "GenbankFile", which can be
instantiated with either file_contents or file_name as strings. It will
parse the target genbank file and import most or all of the file contents
as properties on the GenbankFile object.
References and Features are imported as lists of dedicated container classes.
The Reference class, GBReference, is pretty dumb, but can attempt to present
the referred-to sequence if needed.
The GBFeature class is where most of the magic happens; this class will store
the "/" delimited metadata from the feature table in a dictionary at self.meta,
and has more advanced code for dealing with the target sequence. It can parse
sequence referral lines like "complement(order(10..40,500..743))" and fetch/
invert/complement the target sequence into a contiguous string accessible by
the property self.sequence.
'''

import nucutils
import re

class GBFeature:
    '''Parser/Container class for genbank feature entries.
    Should be passed a feature block as a list of strings corresponding to lines
    from the feature table, including the "CDS"/"gene"/etc opening word, and a
    GenbankFile object containing the referred-to sequence.
    This object doesn't initially contain the referred-to sequence, but
    rather generates it when its "sequence" property is called, and then
    caches it within the object locally.
    This object is compatible with extended IUPAC, and can attempt to give
    reverse complement for any IUPAC-compatible DNA/RNA sequence, but will
    raise an exception if asked to generate the complement of a hybrid between
    RNA and DNA.'''
    # Genbank file format reference: http://www.insdc.org/files/feature_table.html
    # Includes all the horrible dark magic of feature range specifiers.
    gb_func_finder = re.compile(r'(join|complement|order)(\([0-9<>.,]+\))')
    # using gb_func_finder.findall(string) on the range expression returns
    # output like:
    # [('join', '(12..78,134..202)'), ('complement', '(34..126)')]
    # ..but only functions surrounding NUMBERS are returned, so this only
    # resolves the innermost layer of complex expressions like:
    # complement(join(10..20,100..120))
    # This is be used to render this into:
    # complement(GCATTCTAGGCTAGCTATG), which can be parsed further with:
    gb_func_finder_2 = re.compile(r'(join|complement|order)(\([GATCU]{1,}\))')

    def __init__(self, list_of_lines, parent_genbank_object):
        'Should be provided with the feature block split into a list of strings.'
        self._parent_genbank_object = parent_genbank_object
        self.meta = {}
        self._sequence = ''
        # This is used as a mapping dict to call the right function on
        # complex range specifiers when assembling the sequence:
        self.gb_funcs = {"join":self._gb_join,
                         "complement":self._gb_complement,
                         "order":self._gb_order}
        # This is used as a mapping dict to call the right function on
        # partially-finished hybrid statements where a function call
        # now contains a string of DNA rather than the original range:
        self.hybrid_funcs = {"join":self._join,
                             "complement":self._complement,
                             "order":self._order}

        # Get feature type:
        firstline_bits = list_of_lines[0].strip().split()
        self.type = firstline_bits[0]

        # Replace topline without feature type.
        list_of_lines[0] = ' '.join(firstline_bits[1:])
        # Strip all extra whitespace from every line:
        for line in list_of_lines[:]:
            lineindex = list_of_lines.index(line)
            # We need a unique delimiter to replace "/" in meta area,
            # so first strip the leading "/" in lines beginning with it,
            # then add a more unique "/-meta-/" line.
            line = line.strip()
            if line[0] == "/":
                line = line.lstrip("/")
                line = "/-meta-/" + line
            list_of_lines[lineindex] = line

        # Combine all lines so we can split by delimiter into range:meta1:meta2...
        # Above example becomes: 10..567/gene="ubc42"/number=1
        list_of_lines = ' '.join(list_of_lines)
        list_of_lines = list_of_lines.split("/-meta-/")
        # Above example becomes: ['10..567','gene="ubc42"','number=1']

        span_line = list_of_lines[0]
        self.set_span(span_line)
        for metaval in list_of_lines[1:]:
            # Handles storage of "foo=bar" pairs and exception-catching
            # for oddly formatted things.
            self.store_meta(metaval)

    def set_span(self, spanline):
        '''This stores the sometimes-simple, sometimes-nightmarish sequence range specifier line.
        This only does as much as necessary (usually very little) to prepare the span-line for
        later parsing if the feature's sequence is actually asked for.'''
        # Examples of span-lines:
        # (Previously multiline spanlines will be joined by a space)
        # complement(order(4286215..4286223,4286290..4286292, 4286296..4286301,4286596..4286607,4286611..4286613))
        # complement(4289460..4289591)
        # 4278837..4279160
        spanline = spanline.strip().replace(" ","")
        if "<" in spanline or ">" in spanline:
            spanline = spanline.replace("<","").replace(">","")
            self.fuzzyboundary = True
        self.spanline = spanline

    def store_meta(self, meta_block):
        'Parses a "foo=bar" meta-line. Completed lines are added to self.meta.'
        # TODO: Many meta-keys like db_xref may occur many times per feature!
        # Need to refactor to preserve multiple references, perhaps by converting
        # all metablock entries to lists and appending?
        try:
            # Expected format after above processing is 'translation="HSAGHTCNHAT..."'
            meta_tag_bits = meta_block.split("=")
            meta_name = meta_tag_bits[0].strip().strip(r'/\"()')
            # In case content contained "=" symbol and got split, recombine
            meta_tag_bits = [x.strip().strip(r'"\/()') for x in meta_tag_bits]
            meta_content = '='.join(meta_tag_bits[1:])
            if meta_name == "translation":
                meta_content = ''.join(meta_content.split())
            self.meta[meta_name] = meta_content
        except:
            print(("Error occurred while trying to store following"
                  " feature table metadata:"), meta_block, "..skipping.")

    @property
    def gene(self):
        if "gene" not in self.meta.keys():
            return "None"
        else:
            return self.meta['gene']

    @staticmethod
    def _nativise(gb_range_expression):
        '''Renders genbank ranges into Python lists.
        So this: (x..y,z..w) is converted to: [[x,y],[z,w]].
        This method expects single parenthetised data only, so nested
        data will probably choke it.'''
        outlist = []
        gb_range_expression = gb_range_expression.strip().strip("()")
        for subexpression in gb_range_expression.split(","):
            # We ignore < or >, as they aren't useful for this
            # parser's intended use.
            subexpression.strip("<>")
            # Create a list item for int-coerced ranges, then append to
            # output list.
            subexpression = [int(x) for x in subexpression.split("..")]
            outlist.append([subexpression[0]-1, subexpression[1]])
        return outlist

    def _join(self, joinrange):
        '''Join partially-resolved genbank statements like (TTCAG,GAAC)'''
        joinrange = joinrange.strip().strip("()").split(",")
        return ''.join(joinrange)

    def _gb_join(self, joinrange):
        'Given "(x..y,z..w,[a..b]...)", return "GCGATCTGGCATGCT..."'
        # Convert string-formatted genbank range to list of lists:
        # output is of form [[x,y],[z,w],[a,b]...]
        joinrange = self._nativise(joinrange)
        output_list = []
        for subrange in joinrange:
            # subrange = [int, int]
            # Make a shorter pointer to make code not-ugly:
            search_seq = self._parent_genbank_object.sequence
            output_list.append(search_seq[subrange[0]:subrange[1]])
        return ''.join(output_list)

    def _complement(self, c_string):
        'Given "(GGCCTTCC)", return "GGAAGGCC".'
        # No need to split at commas; genbank format should use join() for
        # compound complements like complement(join(1..3,5..7))
        c_string = c_string.strip().strip("()")
        return nucutils.get_complement(c_string)

    def _order(self, o_string):
        raise NotImplementedError(("Error: Genbank files shouldn't have"
                    " any reason to include order as an outer function.."))

    def _gb_complement(self, c_range):
        'Given (x..y), return "GCAGTTACGTGTACTTTACTGGTG"'
        # Convert string-formatted genbank range to list of lists:
        c_range = self._nativise(c_range)[0]
        # Fetch specified sequence range:
        try: output_sequence = self._parent_genbank_object.sequence[c_range[0]:c_range[1]]
        except IndexError as e:
            print("IndexError:",e,"c_range:",c_range,sep="\n")
        # Get reverse complement:
        output_sequence = nucutils.get_complement(output_sequence)
        return output_sequence

    def _gb_order(self, orderrange):
        'Given (a..b,c..d,e..f), return the range a..f.'
        orderrange = self._nativise(orderrange)
        # Just get the start and finish of the range!
        startbase = orderrange[0][0]
        finishbase = orderrange[len(orderrange)-1][1]
        # Sanity check: is start a smaller number than finish?
        if startbase > finishbase:
            # What to do now? Should we return the complement?
            raise NotImplementedError(("Was asked to fetch an ordered"
                   " range where the start base was larger than the"
                   " finish base: not sure what to do about this!"))
        # Now return the sequence spanned by this range..
        return self._parent_genbank_object.sequence[startbase:finishbase]

    def parse_gb_functions(self, gb_expression):
        'Do first-round parsing of gb functions like complement(join(1..2,5..6))'
        # Find first-round functions to address:
        first_round_functions = self.gb_func_finder.findall(gb_expression)
        for func_call in first_round_functions:
            # Format will be ('', (a..b,c..d))
            func_to_call = self.gb_funcs[func_call[0]]
            # Replace the original expression with the output.
            gb_expression = gb_expression.replace(''.join(func_call[0]),func_to_call(func_call[1]))
        return gb_expression

    @property
    def sequence(self):
        if self._sequence:
            # If the below has been called before, the output will be saved
            # to self._sequence to spare us the trouble next time:
            return self._sequence
        if self.spanline.isnumeric():
            return self._parent_genbank_object.sequence[int(self.spanline)-1]
        returnseq = self.spanline
        # First, check if there's a colon in the span; if so, it's
        # referring to another sequence accession, which we can't handle
        if ":" in returnseq:
            raise NotImplementedError(("This sequence contains a ':'"
                     " character, indicating a reference to another"
                     " sequence accession. This isn't yet implemented."))
        returnseq = returnseq.replace("<","").replace(">","").replace("^","..")
        # self.spanline contains the sequence location specifier from the GB file.
        # Parse for order/complement/join blocks using regex
        first_order_funcs = self.gb_func_finder.findall(self.spanline)
        # Resolve first-order functions, replacing them with the output.
        # Should do nothing if no functions are found.
        if first_order_funcs:
            for func_call in first_order_funcs:
                # func_call is of tuple/str format ("join", "(11..45,60..89)")
                func_resolve_method = self.gb_funcs[func_call[0]]
                # Resolution methods are written to deal with flanking brackets:
                replace_content = func_resolve_method(func_call[1])
                content_to_replace = func_call[0]+func_call[1]
                returnseq = returnseq.replace(content_to_replace, replace_content)
            # Parse output for second-order functions using regex.
            second_order_funcs = self.gb_func_finder_2.findall(returnseq)
            # Resolve second-order functions, replacing them with the output.
            # Should do nothing if no second-order calls were found:
            for func_call in second_order_funcs:
                func_resolve_method = self.hybrid_funcs[func_call[0]]
                replace_content = func_resolve_method(func_call[1])
                content_to_replace = func_call[0]+func_call[1]
                returnseq = returnseq.replace(content_to_replace,replace_content)
        else:
            seq_indices = returnseq.split('.')
            for x in seq_indices[:]:
                if not x: seq_indices.pop(seq_indices.index(x))
            seq_indices = [int(x) for x in seq_indices]
            # Adjust only first index so string indexing works.
            seq_indices[0] = seq_indices[0]-1
            if seq_indices:
                # Resolve the range ['x','y']
                try: returnseq = self._parent_genbank_object.sequence[seq_indices[0]:seq_indices[1]]
                except TypeError as e:
                    print("TypeError with indices: "+str(seq_indices)+" types: [{0}, {1}]".format(type(seq_indices[0]),type(seq_indices[1])),e, sep="\n")
        # Check sequence for any remaining non-nucleotide clutter
        charset = nucutils.deduce_alphabet(returnseq)
        for char in charset:
            if char not in nucutils.iupac_characters:
                raise ValueError("Non-IUPAC character(s) found:"+str(charset))
        # Cache and return completed sequence.
        if self._parent_genbank_object.cache_sequences: self._sequence = returnseq
        return returnseq

class GBReference:
    def __init__(self, list_of_lines, parent_genbank_object):
        'Should be provided with a feature block split into a list of strings.'
        #REFERENCE   6  (bases 1 to 2706)
        #  AUTHORS   Piletz,J.E., Deleersnijder,W., Roth,B.L., Ernsberger,P., Zhu,H. and
        #            Ziegler,D.
        #  TITLE     IRAS splice variants
        #  JOURNAL   Ann. N. Y. Acad. Sci. 1009, 419-426 (2003)
        #   PUBMED   15028621
        #  REMARK    GeneRIF: Results describe three alternatively spliced transcripts
        #            of the human I(1)-imidazoline receptor candidate gene, IRAS.
        self.processors = {"REFERENCE":self.process_topline,
                           "AUTHORS":self.process_authors,
                           "CONSRTM":self.process_consortium,
                           "TITLE":self.process_title,
                           "JOURNAL":self.process_journal,
                           "PUBMED":self.process_pubmed,
                           "REMARK":self.process_remark}
        self._parent_genbank_object = parent_genbank_object
        # Reference number in the genbank file:
        self.refno = None
        # Range of the sequence that this reference applies to:
        self.range = None
        self.seqrange = None
        self.authors = 'None'
        self.consortium = 'None'
        self.title = 'None'
        self.journal = 'None'
        self.pubmed = 'None'
        self.remark = 'None'
        subsection = []
        for line in list_of_lines:
            firstword = line.strip().split()[0]
            if firstword in self.processors.keys():
                # This new line is a header, so prior header (if any) is finished.
                if subsection:
                    # If there is a prior section in the buffer, process
                    # it and empty the buffer. Get first word of buffer:
                    prior_firstword = subsection[0].strip().split()[0]
                    # Use first word to call the appropriate processor method:
                    self.processors[prior_firstword](' '.join(x.strip() for x in subsection).strip())
                    # Clear the buffer:
                    subsection = []
            # Append new line to the buffer.
            subsection.append(line)
        # After the for-loop the last block will still be in the buffer, so
        # process that now.
        firstword = subsection[0].strip().split()[0]
        self.processors[firstword](' '.join(x.strip() for x in subsection).strip())

    def process_topline(self, line):
        #REFERENCE   6  (bases 1 to 2706)
        # Not all reference toplines contain (bases x to y):
        if " (bases" in line:
            # Cut off the closing bracket if any:
            line = line.strip(")")
            line_bits = line.split("(")
            # line should now be REFERENCE   X, same as if (bases x to y) weren't
            # there in the first place; makes code below this if-block easier.
            line = line_bits[0].strip()
            # Strip all the characters we expect from the edges ("bases ")
            span_line = line_bits[1].strip("abes ")
            span_line = span_line.split(" to ")
            if " " in span_line[1]:
                # If there's an unprocessed line after the range, it can
                # lead to a crash when trailing stuff gets coerced to int()
                # below.. so make sure any trailing crap gets sliced off.
                extra_crap = span_line[1].split()
                span_line[1] = extra_crap[0]
                print("Unprocessed data found after range:",' '.join(extra_crap))
            self.range = (int(span_line[0]), int(span_line[1]))
            # self.range is the numbers as a human might read them, but
            # self.seqrange is used for string indexing of the actual sequence.
            self.seqrange = [x-1 for x in self.range]
        # Either way, line now just contains "REFERENCE   XX"..
        line = line.split()
        self.refno = int(line[1])

    @property
    def sequence(self):
        'This property returns the sequence referred to (if any) in the REFERENCE topline.'
        if self.seqrange:
            return self._parent_genbank_object['Sequence'][self.seqrange[0]:self.seqrange[1]]

    def process_authors(self, line):
        if line[:7] == "AUTHORS":
            line = line[7:].lstrip()
        self.authors = line

    def process_consortium(self, line):
        if line[:7] == "CONSRTM":
            line = line[7:].lstrip()
        self.consortium = line

    def process_title(self, line):
        if line[:5] == "TITLE":
            line = line[5:].lstrip()
        self.title = line

    def process_journal(self, line):
        if line[:7] == "JOURNAL":
            line = line[7:].lstrip()
        self.journal = line

    def process_pubmed(self, line):
        if line[:6] == "PUBMED":
            line = line[6:].lstrip()
        self.pubmed = line

    def process_remark(self, line):
        if line[:6] == "REMARK":
            line = line[6:].lstrip()
        self.remark = line

class GenbankFile(dict):
    '''Parses a genbank-formatted string or file.
    Must be provided either of the "file_contents" or "file_name" arguments.
    file_name should be a path to the target file,
    file_contents should be a string as if file.read() directly from a file.'''
    def __init__(self, file_contents=None, file_name=None, cache=True):
        'Accepts either a genbank filename or contents of same.'
        self.block_parsers = {"ORIGIN":self.process_sequence,
                              "FEATURES":self.process_features,
                              "COMMENT":self.process_comment,
                              "VERSION":self.process_version,
                              "ACCESSION":self.process_accession,
                              "DEFINITION":self.process_definition,
                              "REFERENCE":self.process_reference,
                              "SOURCE":self.process_source,
                              "DBLINK":self.process_dblink,
                              "KEYWORDS":self.process_keywords,
                              "LOCUS":self.process_locus,
                              "//":self._ignore}
        self['Sequence'] = ''
        self['References'] = []
        self['Accession'] = []
        self['Definition'] = ''
        self['Locus'] = ''
        self['Source'] = ''
        self['Metadata'] = {"Keywords":[],
                            "DBLink":{},
                            "Version":'',
                            "Comment":'',}
        self["Features"] = []
        if not file_contents:
            if not file_name:
                raise ValueError("Either a genbank filename or a string with genbank-formatted data must be provided.")
            if not isinstance(file_name, str):
                raise ValueError("Filename must be provided as a string.")
            try:
                with open(file_name) as Genbank_File:
                    file_contents = Genbank_File.read()
            except IOError as e:
                print("An exception occurred while trying to open file '{0}':".format(file_name), e)
            self.extract_indent_blocks(file_contents)
            self.process_indent_blocks()

        # Directs subordinate GBFeature objects whether or not they should
        # retain a copy of their parsed/converted sequences in memory. If
        # set to True, then subsequent lookups will be faster, but at cost
        # of RAM. Some uses of GenbankFile may call sequence for every object
        # which will incur a processing cost up-front anyway, and in these
        # cases it is a matter of taste whether to preserve future processing
        # power or RAM.
        self['cache_sequences'] = cache

    def process_version(self, version_strings):
        'Saves version to self.metadata["Version"].'
        version_strings[0] = ' '.join(version_strings[0].split()[1:])
        version_strings = [x.strip() for x in version_strings]
        self['Metadata']['Version'] = ' '.join(version_strings)

    def process_comment(self, comment_strings):
        'Saves comment to self.metadata["Comment"].'
        comment_strings[0] = ' '.join(comment_strings[0].split()[1:])
        comment_strings = [x.strip() for x in comment_strings]
        self['Metadata']['Comment'] = ' '.join(comment_strings)

    @staticmethod
    def _count_indent(line_string):
        'Compares string to lstripped version of itself to determine indent level.'
        return len(line_string) - len(line_string.lstrip())

    def _ignore(self, fooline):
        pass

    def process_locus(self, locus_strings):
        'Currently just adds the locus line to the self["Locus"] key.'
        full_locus_line = ' '.join(locus_strings)
        self['Locus'] = ' '.join(full_locus_line.split()[1:])

    def process_keywords(self, keywords_strings):
        'This splits keywords, if present, and appends to Metadata "Keywords" list.'
        # Empty:     KEYWORDS    .
        # Not Empty: KEYWORDS    complete genome.
        for line in keywords_strings:
            line = line.strip().strip(".")
            if line[:8] == "KEYWORDS":
                if len(line) > 8:
                    line = line[8:]
                else:
                    # Nothing but "KEYWORDS" left after stripping. Ignore.
                    break
            # Add all elements of line.split() to the keywords list.
            self['Metadata']['Keywords'].extend(line.split())

    def process_dblink(self, dblink_strings):
        'For each database link, this adds a dictionary key for the database with an accession value.'
        # DBLINK      BioProject: PRJNA78
        #             Sequence Read Archive: ERS000041
        for line in dblink_strings:
            line = line.strip()
            if "DBLINK" in line:
                line = line[6:].strip()
            linebits = line.split(":")
            self['Metadata']['DBLink'][linebits[0].strip()] = linebits[1].strip()

    def process_definition(self, definition_strings):
        # DEFINITION  Homo sapiens nischarin (NISCH), transcript variant 3, mRNA.
        self['Definition'] = ' '.join(definition_strings[0].strip().split()[1:])

    def process_accession(self, accession_strings):
        # Simple Format: ACCESSION   NM_001276294
        # Annoying Format: ACCESSION   AE016825 AE016910 AE016911 AE016912 AE016913 AE016914 AE016915
        #                              AE016916 AE016917 AE016918 AE016919 AE016920 AE016921 AE016922
        #                              AE016923 AE016924 AE016925
        # Join all substrings by a space character:
        accessions_list = ' '.join(accession_strings)
        # Use 2nd word onwards as accessions list:
        self['Accession'] = accessions_list.split()[1:]

    def process_source(self, source_strings):
        # Format: SOURCE      Homo sapiens (human)
        # Note my hideous method of removing the first word. Yes, cringe, it won't help.
        self['Source'] = ' '.join(source_strings[0].strip().split()[1:])

    def extract_indent_blocks(self, gbfile_contents):
        'Parses a genbank file and extracts indented blocks into class properties.'
        gb_lines = gbfile_contents.strip().splitlines()
        self.indent_blocks = []
        this_block = []
        # This for/if construct divides a genbank file into indented blocks:
        for line in gb_lines:
            if not line:
                continue
            if not line.strip():
                # Skip empty lines
                continue
            if line[0] != " " and this_block:
                # Add the list of lines to the list of sub-blocks.
                self.indent_blocks.append(this_block)
                this_block = []
            this_block.append(line)
        if this_block:
            self.indent_blocks.append(this_block)

    def process_indent_blocks(self):
        'Processes each indented block of a Genbank file with a sub-parser.'
        for block in self.indent_blocks:
            # Remember: each "block" is a *list* of lines from that block.
            # Detect first word of first block line
            firstword = block[0].split()[0].strip().upper()
            # Use first word as key to select appropriate sub-parser,
            # then pass the sub-block to sub-parser.
            try:
                self.block_parsers[firstword](block)
            except KeyError: # No parser defined for this block
                print("No parser found for block '{0}'.".format(firstword))
                continue
            # Done! Sub-parsers are written to add output to self.data_blocks.

    def process_sequence(self, origin_block_lines):
        'Expects the sequence from the end of a Genbank file, from "ORIGIN" onwards.'
        # Rather than doing string contatenation (i.e. 'this' =+ 'that'),
        # we build a list of substrings and then ''.join(list) them, as this
        # is far faster and memory efficient.
        output = []
        for line in origin_block_lines:
            if not line.strip():
                # Ignore empty lines.
                continue
            elif line.strip() == "//":
                # Stop processing at end of file indicator.
                break
            elif line.strip() == "ORIGIN":
                # Skip the "ORIGIN" line.
                continue
            else:
                # Append a line to the output list, minus any nucleotide-numbers or whitespace.
                line_content = line.strip().strip("1234567890").replace(" ", "")
                output.append(line_content.upper())
        self['Sequence'] = ''.join(output)

    def process_features(self, features_block_lines):
        'Parses through features block to extract genes, CDS, mRNA etc.'
        def addfeature(new_feature):
            NewFeature = GBFeature(this_feature, self)
            self['Features'].append(NewFeature)
        this_feature = []
        for line in features_block_lines[1:]:
            indent_lvl = self._count_indent(line)
            if indent_lvl == 5:
                # Subsections should be at this indent
                if this_feature:
                    addfeature(this_feature)
                    this_feature = []
            this_feature.append(line)
        if this_feature:
            addfeature(this_feature)

    def process_reference(self, reference_strings):
        'Parse a reference block and append to list.'
        try:
            self['References'].append(GBReference(reference_strings, self))
        except:
            print("An error occurred while processing reference:\n\t","\n\t".join(reference_strings))

    def __getattr__(self, attribute):
        if attribute in self.keys():
            return self[attribute]
        elif attribute.title() in self.keys():
            # Check if the string is in self.keys(), but in title-case:
            return self[attribute.title()]
        else:
            raise AttributeError("Attribute {0} neither in object namespace nor in genbank object dictionary.".format(attribute))

def testfeatures(gb_obj):
    errors = []
    feature_types = []
    meta_keys = []
    excessive_meta_features = []
    for feature in gb_obj.features:
        try:
            x = feature.sequence
        except:
            errors.append(feature)
            print("Error found in feature:",feature.type,feature.spanline)
        if feature.type not in feature_types:
            feature_types.append(feature.type)
        for key in feature.meta.keys():
            if key not in meta_keys:
                meta_keys.append(key)
            if len(key) > 15:
                excessive_meta_features.append(feature)
    return errors, feature_types, meta_keys, excessive_meta_features

def tests():
    print("Importing test genomes, may take several seconds..")
    testgenomes = []
    testgenomes.append(GenbankFile(file_name="test_genomes/E.coli_DH10B.gb"))
    testgenomes.append(GenbankFile(file_name="test_genomes/E.coli_K12_MG1655.gbk"))
    testgenomes.append(GenbankFile(file_name="test_genomes/E.coli_K12_W3110.gbk"))
    testgenomes.append(GenbankFile(file_name="test_genomes/C.difficile_630_uid78.gb"))
    testgenomes.append(GenbankFile(file_name="test_genomes/Cyanobacteria_bacterium_Yellowstone_A-Prime_uid16251.gb"))
    print("Imported test genomes:\n\t", "\n\t".join([x['Definition'] for x in testgenomes]))
    def uniquify(sequence):
        'Raplidly reduces a sequence to unique-components only, in order.'
        seen = set()
        return [x for x in sequence if x not in seen and not seen.add(x)]
    overall_features = []
    overall_meta = []
    for genome in testgenomes:
        print("Testing sequence:",genome.source)
        errors, feature_types, meta_keys, big_meta = testfeatures(genome)
        print("\tEncountered feature types:","; ".join(feature_types))
        print("\tEncountered meta keys:", "; ".join(meta_keys))
        overall_features.extend(feature_types)
        overall_meta.extend(meta_keys)
    print("\nMeta encountered across test genomes:\n",'; '.join(uniquify(overall_meta)))
    print("\nFeature types encountered across test genomes:\n","; ".join(uniquify(overall_features)))

if __name__ == "__main__":
    tests()

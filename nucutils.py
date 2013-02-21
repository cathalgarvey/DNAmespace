'''nucutils - Nucleotide-wrangling functions and constants.
by Cathal Garvey
Part of the DNAmespace project. License accessible as nucutils.license.
'''
from gnulicenses import Affero as license

iupac_characters = [ 'A', 'T', 'C', 'G', 'U',    # Canonical bases
                     'B', 'V', 'D', 'H',     # B=Not A, V=Not T, D=Not C, H=Not G
                     'S', 'W',               # Strong and Weak: GC vs. AT
                     'K', 'M',               # "Keto" and "aMino": GT vs. AC
                     'R', 'Y',               # "puRine" and "pYrimidine": AG vs CT
                     'N', '.', '-' ]         # Wildcards: any N.

# All complements are in lowercase so that string substitution is simple:
# just set string to uppercase (if not already), replace all characters
# with lowercase complements, then if desired set resulting string to
# uppercase again.
dnacomplement = {"A":"t","T":"a","G":"c","C":"g"}
rnacomplement = {"A":"u","U":"a","G":"c","C":"g"}
dnaiupaccomplement = {   # Translates to an IUPAC sequence of potential complements.
              "A":   "t", "T": "a",
              "C":   "g", "G": "c",
              "W":   "w", "S": "s",
              "B":   "v", "V": "b",
              "H":   "d", "D": "h",
              "M":   "k", "K": "m",
              "R":   "y", "Y": "r",
              "V":   "b", "B": "v",
              "N":   "n"}
rnaiupaccomplement = {   # Translates to an IUPAC sequence of potential complements.
              "A":   "u", "U": "a", "C":   "g", "G": "c",
              "W":   "w", "S": "s", "B":   "v", "V": "b",
              "H":   "d", "D": "h", "M":   "k", "K": "m",
              "R":   "y", "Y": "r", "V":   "b", "B": "v",
              "N":   "n"}

def _uniquify(string):
    '''Reduces a string down to its component characters.
    This is a fast, order-preserving function for removing duplicates
    from a sequence/list, by "Dave Kirby"
    Found here: http://www.peterbe.com/plog/uniqifiers-benchmark'''
    seen = set()
    return [x for x in string if x not in seen and not seen.add(x)]

def deduce_alphabet(string):
    charset = _uniquify(string)
    for character in charset:
        if character not in iupac_characters:
            raise ValueError(("Non-base found in nucleotides string."
                  " The string consists of the following characters:\n")
                  +str(charset))
    if "T" in charset and "U" in charset:
        raise ValueError(("Cannot get reverse complement of a hybrid"
                          " DNA/RNA sequence."))
    elif "U" in charset:
        basedict = rnaiupaccomplement
    else:
        basedict = dnaiupaccomplement
    return basedict

def get_complement(nucleotides):
    'Given a string of nucleotides (RNA *or* DNA), return reverse complement.'
    # Reverse the string:
    nucleotides = nucleotides[::-1]
    # Determine molecule type:
    basedict = deduce_alphabet(nucleotides)
    for base in basedict.keys():
        # Replace each base with its complementary base in lowercase:
        # (lowercase means previously substituted bases aren't replaced)
        nucleotides = nucleotides.replace(base, basedict[base])
    # Return sequence in uppercase:
    return nucleotides.upper()

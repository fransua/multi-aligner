#!/usr/bin/python
"""
30 Sep 2011


"""

__author__  = "Francois-Jose Serra"
__email__   = "francois@barrabin.org"
__licence__ = "GPLv3"
__version__ = "0.0"



from itertools               import groupby
from align.utils.seq_utils import translate
from sys                     import stderr

def parse_fasta(path, genetic_code=None, store='seq'):
    """
    """
    def check_header(line):
        return line.startswith('>')
    def get_sep (line):
        if ' |' in line:
            return ' |'
        elif '\t|' in line:
            return '\t|'
        elif '|' in line:
            return '|'
        elif '\t' in line:
            return '\t'
        elif ' ' in line:
            return ' '
        else:
            return '\n'
    sep = None
    seqs = {}
    with open(path) as fasta:
        for is_header, lines in groupby (fasta, key=check_header):
            text = ''.join (line.strip() for line in lines)
            if is_header:
                if not sep:
                    sep = get_sep(text)
                if sep == '\n':
                    head, descr = text[1:], ''
                else:
                    head, descr = text[1:].split(sep, 1)
                seqs [head] = {'descr': descr}
            else:
                seq = text.replace('\n', '')
                if genetic_code:
                    if len (text)%3 != 0:
                        seq = seq[:-(len(seq)%3)]
                    try:
                        seqs[head]['prot']  = translate(text, genetic_code,
                                                       stop=True)
                        seqs[head]['codon'] = [seq[i:i+3] for i in xrange(0,len(seq), 3)]
                    except KeyError: # in case sequence == 'Sequence unavailable'
                        del (seqs [head])
                        print >> stderr, 'No sequence found for ' + head
                        print >> stderr, text
                        continue
                seqs [head][store] = seq
    return seqs



def parse_mcoffee_score(path, sequences):
    """
    """
    def is_header(line):
        return line=='\n'
    with open (path) as score:
        for is_header, lines in groupby (score, key=is_header):
            if is_header: continue
            for line in lines:
                if line.startswith ('T-COFFEE,'): break
                if line.startswith ('cons  '): continue
                name, sco = line.split ()
                sco = ['0' if n.isalpha() else n for n in sco]
                sequences[name].setdefault ('score', []).extend (sco)


def parse_mcoffee_aln(path, sequences):
    """
    """
    def is_header(line):
        return line=='\n'
    with open (path) as score:
        for is_header, lines in groupby (score, key=is_header):
            if is_header: continue
            for line in lines:
                if line.startswith ('CLUSTAL'): break
                if line.startswith ('    '): continue
                name, seq = line.split ()
                sequences[name].setdefault ('aa_ali', []).extend (seq)

#!/usr/bin/python
"""
11 Jul 2011


"""

__author__  = "Francois-Jose Serra"
__email__   = "francois@barrabin.org"
__licence__ = "GPLv3"
__version__ = "0.0"

from itertools import groupby
from re import sub
from re import compile as compil
from sys import stdout
from align.utils.seq_utils import translate, get_genetic_code

class Sequences():
    def __init__(self, genetic_code=None, path=None):
        self.genetic_code = get_genetic_code(genetic_code)
        self.__sequences = {}
        self.headers = []
        self.items = []
        if path:
            self.read(path)
        
    def __getitem__(self, name):
        return self.__sequences[name]

    def __setitem__(self, name, value):
        self.__sequences[name] = value

    def __delitem__(self, name):
        del (self.__sequences[name])
        self.headers.pop(self.headers.index(name))

    def __iter__(self):
        for seq in self.__sequences:
            yield seq

    def __repr__(self):
        seqs = ''
        for head in self.headers:
            seqs += '>%s |%s\n' % (head, self[head]['descr'])
            for item in self.items:
                seqs += '    %s length: %s\n' % (item, len(self[head][item]))
        return seqs

    def __str__(self):
        seqs = ''
        for head in self.headers:
            seqs += '>%s |%s\n' % (head, self[head]['descr'])
            for item in self.items:
                seqs += '%s: %s\n' % (item, self[head][item])
        return seqs

    def write(self, outfile=None, item='seq', reverse=False, width=60, descr=False):
        """
        Write sequence object to file in fasta format
        
        :argument None outfile: path to outfile, if None than, print to stdout
        :argument seq item: what to put in place of sequence
        :argument False reverse: wether to reverse or not the sequence
        :argument 60 width: number of sites per line when printing sequence
        :argument False descr: put description of sequence also, not recommended if you are not sure how the aligner will read it.
        
        """
        if outfile:
            out = open (outfile, 'w')
        else:
            out = stdout
        wsub = compil('([A-Za-z-]{'+str(width)+'})')
        for elt in self:
            if decr:
                out.write ('>%s |%s\n' % (elt, self[elt]['descr']))
            else:
                out.write ('>%s\n' % (elt))
            seq = self[elt][item][::-1] if reverse else self[elt][item]
            seq = seq if type(seq) is str else ''.join(seq)
            out.write ('%s\n' % (sub(wsub, '\\1\n', seq)))
        if outfile:
            out.close()
        
    
    def read(self, path, genetic_code=None, store='seq'):
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
                    self [head] = {'descr': descr}
                else:
                    seq = text.replace('\n', '')
                    if genetic_code:
                        if len (text)%3 != 0:
                            seq = seq[:-(len(seq)%3)]
                        try:
                            self[head]['prot']  = translate(text, genetic_code,
                                                           stop=True)
                            self[head]['codon'] = [seq[i:i+3] for i in xrange(0,len(seq), 3)]
                        except KeyError: # in case sequence == 'Sequence unavailable'
                            del (self[head])
                            print >> stderr, 'No sequence found for ' + head
                            print >> stderr, text
                            continue
                    self [head][store] = seq
                    self.headers.append(head)
        for item in self[head]:
            if item == 'descr': continue
            self.items.append(item)
        
    
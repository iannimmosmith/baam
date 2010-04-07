''' Testing BANARM module '''

import numpy as np

from nose.tools import assert_true, assert_false, assert_equal

from numpy.testing import assert_array_equal, assert_array_almost_equal

import BANARM as baam

def test_map_lists():
    def f(x):
        return [2*x+1,2*x]
    yield assert_equal, baam.map_lists(f,[0,3,7]), [15, 14, 1, 0, 7, 6]
    
def test_uniquify():
    yield assert_equal, baam.uniquify([('mu', 'riyH', 'PFX[m] '), ('mu', 'riyH', 'PFX[m] V/3 '), ('mu', 'riyH', 'PFX[m] V/3 ')]), [('mu', 'riyH', 'PFX[m] '), ('mu', 'riyH', 'PFX[m] V/3 ')]

def test_analyst():
    yield assert_equal, baam.analyst(('mu', 'riyH', 'PFX[m] ')), [('mu', 'ryH', 'fiEl', 'PFX[m] C/3 3R ')]
    yield assert_equal, baam.analyst(('','DAfiy','')), [('', 'Dfy', 'fAEil', 'C/5 3R '), ('', 'DAy', 'fEfil', 'C/5 3R '), ('', 'DfY', 'fAEil', 'C/5 3R RootEnd[wyY] '), ('', 'DAfw', 'fElil', 'C/5 4R RootEnd[wyY] '), ('', 'DAf', 'fEliy', 'C/5 3R '), ('', 'DAfY', 'fElil', 'C/5 4R RootEnd[wyY] '), ('', 'Afy', 'DfEil', 'C/5 3R '), ('', 'DAY', 'fEfil', 'C/5 3R RootEnd[wyY] '), ('', 'DAf', 'fEliy', 'V/5 3R '), ('', 'DAfy', 'fElil', 'C/5 4R '), ('', 'AfY', 'DfEil', 'C/5 3R RootEnd[wyY] '), ('', 'DAw', 'fEfil', 'C/5 3R RootEnd[wyY] '), ('', 'Afw', 'DfEil', 'C/5 3R RootEnd[wyY] '), ('', 'Dfw', 'fAEil', 'C/5 3R RootEnd[wyY] ')]

def test_parse_list():
    yield assert_equal, baam.parse_list([('pa','sa','ca'),('pb','sb', 'cb')]), [('pa','sa','ca'),('pb','sb','cb')] 
    
def test_parse_prefix():
    yield assert_equal, baam.parse_prefix(('muriyH','')), [('', 'muriyH', ''), ('mu', 'riyH', 'PFX[m] ')]

                       
def test_substitute_initial():
    yield assert_equal, baam.substitute_initial([('','>us~','')]), [('', '>us~', ''), ('', 'Aus~', 'IN '), ('', "'us~", 'IN '), ('', '<us~', 'IN '), ('', '{us~', 'IN ')]

def test_geminate_tildes():
    yield assert_equal, baam.geminate_tildes([('','>us~','')]),[('', '>us~', ''), ('', '>uss', 'GM/3 ')]


def test_find_tilde_locations():
    yield assert_equal, baam.find_tilde_locations('abc~de~'), [3, 6]

def test_powerset_graycode():
    yield assert_equal, [list(p) for p in baam.powerset_graycode([3,6])], [[], [3], [3, 6], [6]]

def test_two_consonant_expand():
    yield assert_equal, baam.two_consonant_expand([('','>us~','')]), [('', '>us~', ''), ('', '>ius~', 'CC[i] '), ('', '>yus~', 'CC[y] '), ('', '>uss~', 'CC2 ')]
    yield assert_equal, baam.two_consonant_expand([('','>us~',''), ('','>ud~','')]), [('', '>us~', ''), ('', '>ius~', 'CC[i] '), ('', '>yus~', 'CC[y] '), ('', '>uss~', 'CC2 '), ('', '>ud~', ''), ('', '>iud~', 'CC[i] '), ('', '>yud~', 'CC[y] '), ('', '>udd~', 'CC2 ')]

def test_meta_expand_list():
    yield assert_equal, baam.meta_expand_list([('','DAfiy','')]), [('', 'D%fi#',  'C/5 '), ('' ,'D%fiy', 'V/5 ')]

def test_meta_expand():
    yield assert_equal, baam.meta_expand(('','DAfiy','')), [('', 'D%fi#',  'C/5 '), ('' ,'D%fiy', 'V/5 ')]

'''
def test_dictionary_sets():
    #Create parsed dictionary
    parsed_dictionary = baam.make_parse_dictionary()
    #Create frequency dictionaries for stems, roots and patterns
    stem_dictionary = baam.make_BANARM_dictionary('stem_pointed_token')
    root_dictionary = baam.make_BANARM_dictionary('root_token')
    pattern_dictionary = baam.make_BANARM_dictionary('pattern_token')
    answer = baam.dictionary_sets(stem_dictionary, parsed_dictionary)
    yield assert_equal, answer['commentary'], 'There are 33651 stems in the dictionary\nThere are 38270 stems in the corpus\nThere are 18795 stems unique to the dictionary\nThere are 23414 stems unique to the corpus\nThere are 14856 stems common to both corpus and dictionary\n' 

    

'''

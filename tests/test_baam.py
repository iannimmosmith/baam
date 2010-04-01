''' Testing BANARM module '''

import numpy as np

from nose.tools import assert_true, assert_false, assert_equal

from numpy.testing import assert_array_equal, assert_array_almost_equal

import BANARM as baam


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
    yield assert_equal, baam.two_consonant_expand([('','>us~','')]), [[('', '>us~', '')], [('', '>ius~', 'CC[i] '), ('', '>yus~', 'CC[y] '), ('', '>uss~', 'CC2 ')]]

def test_meta_expand():
    yield assert_equal, baam.meta_expand([('','DAfiy','')]), set([('', 'D%fiy', 'V/5 '), ('' ,'D%fi#', 'C/5 ')])
     

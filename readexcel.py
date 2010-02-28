#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Author: Ian Nimmo-Smith
Description: Wrapper for 'xlrd' Python library to import and read Excel files
'''

try:
  import xlrd
except ImportError:
  print('xlrd is not installed.')

##try:
##    import xlwt
##except ImportError:
##    print('xlwt is not installed.')

##import scipy as sp
##from scipy import linalg as lg
##from scipy import io
##import numpy as np
##import os
##import string
##import struct
##import gzip
##from optparse import OptionParser

def loadexcel(filename):
  '''
  Load the full excel spreadsheet book using xlrd
  book has the following important properties 

  book.nsheets
  book.sheet_names()
  sh = book.sheet_by_index(0)  
  sh = book.sheet_by_name(sheetname)
  
  sh.name, sh.nrows, sh.ncols
  print "Cell D30 is", sh.cell_value(rowx=29, colx=3)
  for rx in range(sh.nrows):
    print sh.row(rx)
  
  More about xlrd 
  http://www.lexicon.net/sjmachin/xlrd.html
  And an example of usage here
  http://www.lexicon.net/sjmachin/README.html
  
  '''
  try:
    return xlrd.open_workbook(filename)
  except:
    print('Cannot open file.')
    pass
    
##def saveexcel(filename):
##    
##    '''
##    Similar to xlrd but writing this time
##    http://pypi.python.org/pypi/xlwt
##    '''
##    try:
##        xlwt.Workbook
##        book=Workbook()
##        sheet1 = book.add_sheet('Sheet 1')
##        sheet1.write(0,0, 'Words')
##        book.save(filename)
##        
##    except:
##        print('Cannot create workbook.')
    

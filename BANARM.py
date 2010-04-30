'''
Author: Ian Nimmo-Smith
Description: Python library for Bayesian ANalysis of ARabic Morphology (BANARM)
Date: 28 Apr 2010
'''

import readexcel as form
import string
import cPickle
import numpy as np
import re

slot_names3=['f', 'E', 'l']
slot_names4=['f', 'E', 'l', 'l']
initials = "A'<>{"
glottal_stops = "'|>&<}`{"
glottal_final_roots = "'>A<{}"
glottal_medial_roots = "&A'>"
extra_letters = "'Atslmnwyo"
pointing = 'aiu~oN'
vowels = 'aeiou~Awy'
tilde = '~'

baampath = r'/home/ian/Sami/baam/'

def parse_list(prefix_stem_code_list):
    return [parse(outer) for outer in prefix_stem_code_list]

def parse(p_s_c):
    return p_s_c    

def gen_comb(p_s_c):
    lis = []
    (prefix,stem,raw_code)=p_s_c
    places=find_nonvowel_locations(stem)
    if len(places) >= 3:
        code = raw_code
        if '3R' not in code:
            code+='3R '
        for root in generate_combinations(places,3):
            try_pattern=list(consonant_meta(stem))
            try_root=string.join([try_pattern[x] for x in root],'')
            for x in range(3):
                try_pattern[root[x]]=slot_names3[x]
            try_pattern=string.join(try_pattern,'')
            lis += [(prefix, try_root, try_pattern, code)]
    if len(places) >= 4:
        code = raw_code
        if '4R' not in code:
            code+='4R '
        for root in generate_combinations(places,4):
            try_pattern=list(consonant_meta(stem))
            try_root=string.join([try_pattern[x] for x in root],'')
            for x in range(4):
                try_pattern[root[x]]=slot_names4[x]
            try_pattern=string.join(try_pattern,'')
            lis += [(prefix, try_root, try_pattern, code)]
    return lis

def generate_combinations(items, n):
    '''
    Combine sets of size n from items
    Perhaps this should go in a combinatorics utility module?
    '''
    if n == 0:
        yield []
    else:
        for i in xrange(len(items)):
            for cc in generate_combinations(items[i+1:], n-1):
                yield [items[i]] + cc

def root_end(prefix_root_pattern_code, compstring, addcode):
    result  = [prefix_root_pattern_code]
    prefix, root, pattern, code =  prefix_root_pattern_code
    if root[-1] in compstring:
        for char in compstring:
            if root[-1] != char:
                result += [(prefix, root[:-1]+char,pattern,code+addcode+' ')]
    return result
                
def pattern_end(prefix_root_pattern_code, compstring, addcode):
    result  = [prefix_root_pattern_code]
    prefix, root, pattern, code =  prefix_root_pattern_code
    if pattern[-1] in compstring:
        for char in compstring:
            if pattern[-1] != char:
                result += [(prefix, root ,pattern[:-1]+char, code+addcode+' ')]
    return result
    
def pattern_start(prefix_root_pattern_code):
    result  = [prefix_root_pattern_code]
    prefix, root, pattern, code =  prefix_root_pattern_code
    if pattern[0] in initials:
        for char  in initials:
            if char != pattern[0]:
                result += [(prefix, root, char+pattern[1:], code+'PattInitGlot ')]
    return result

def glottal_medial_root(prefix_root_pattern_code):
    result  = [prefix_root_pattern_code]
    prefix, root, pattern, code =  prefix_root_pattern_code
    if len(root) == 3 and root[1] in glottal_medial_roots:
        for char in glottal_medial_roots:
            if root[1] != char:
                result += [(prefix,root[0]+char+root[2],pattern, code+'GlotMedRoot ')]
    return result

def glottal_final_root(prefix_root_pattern_code):
    result = [prefix_root_pattern_code]
    prefix, root, pattern, code = prefix_root_pattern_code
    if len(root) == 3 and root[2] in glottal_final_roots:
        for char in glottal_final_roots:
            if root[1] != char:
                result += [(prefix,root[:2]+char,pattern,code+'GlotFinRoot ')]
    return result

def cvc_root(prefix_root_pattern_code):
    prefix, root, pattern, code = prefix_root_pattern_code
    last=pattern.find('E')
    cvc = [prefix_root_pattern_code]
    if root[1] == 'A':
        cvc += [(prefix, root[0]+"y"+root[2], pattern[:last]+"aEa"+pattern[last+1:], \
                     code+'Sub[A/(y,aEa)] ')]
        cvc += [(prefix, root[0]+"w"+root[2], pattern[:last]+"aEa"+pattern[last+1:], \
                     code+'Sub[A/(w,aEa)] ')]
    elif root[1] == 'w':
        cvc += [(prefix, root, pattern[:last]+"uEu"+pattern[last+1:], \
                     code+'Sub[w/(w,uEu)] ')]
    elif root[1] == 'y':
        cvc += [(prefix, root, pattern[:last]+"iEi"+pattern[last+1:], \
                     code+'Sub[y/(y,iEi)] ')]
    return cvc
    

####    oddsandends = set([])
####    for (root,pattern) in parts:
####        if root[2] in "Yyw":
####            for ending in "Yyw":
####                oddsandends.add((root[0:2]+ending,pattern))
####        if root[1] in "Awy":
####            for middle in "Awy":
####                oddsandends.add((root[0]+middle+root[2],pattern))
####        if root[2] in glottal_stops:
####            for ending in glottal_stops:
####                oddsandends.add((root[0:2]+ending,pattern))
####    parts.update(oddsandends)
##    return parts

def analyst(prefix_stem_code):
    psc = two_consonant_expand([prefix_stem_code])
    #psc_m = sum(map(meta_expand, psc_c),[])
    psc = map_lists(geminate_tildes, psc)
    psc = map_lists(meta_expand, psc)
    psc = map_elements(parse,psc)
    prpc = map_lists(gen_comb,psc)
    #prpc1 = map_lists(root_end_wyY,prpc)
    #prpc = map_lists(cvc_root,prpc)
    prpc = map_lists_args(root_end, (prpc,'wyY', 'RootEnd[wyY]'))
    prpc = map_lists_args(pattern_end, (prpc,'wyY','PattEnd[wyY]'))
    prpc = map_lists_args(pattern_end, (prpc, "&'}","PattEnd[&'}]"))
    prpc = map_lists_args(pattern_start, (prpc,))
    prpc = map_lists_args(glottal_medial_root, (prpc,))
    prpc = map_lists_args(glottal_final_root, (prpc,))
    return prpc


def map_lists(fun,lis):
    # fun is a function which returns lists
    # lis is a list

    #function to flatten and uniquify
    #the application of a function to a list
    #It is assumed that fun returns lists
    #return list(set(sum([[apply(fun,[el]) for el in lis]],[])))
    return uniquify(sum([m for m in sum([[apply(fun, [el]) for el in lis]], [])], []))
    #uniquify(sum([[apply(fun,[el]) for el in lis]],[]))

def map_lists_args(fun,argtuple):
    # fun is a function which returns lists
    # argtuple is a tuple with first element a list
    # for fun to be evaluated on with additional arguments
    # from the subsequent elements of argtuple (if any)

    #function to flatten and uniquify
    #the application of a function to a list
    #It is assumed that fun returns lists
    #return list(set(sum([[apply(fun,[el]) for el in lis]],[])))
    lis = argtuple[0]
    otherargs = argtuple[1:]
    lists=sum([m for m in sum([[fun(*((el,)+otherargs)) for el in lis]], [])], [])
    return uniquify(lists)
    #uniquify(sum([[apply(fun,[el]) for el in lis]],[]))

def map_elements(fun,lis):
    # fun is a function which returns values
    # lis is a list

    #function to flatten and uniquify
    #the application of a function to a list
    #It is assumed that fun returns the values we want to uniquely listify
    #return list(set(sum([[apply(fun,[el]) for el in lis]],[])))
    return uniquify(sum([[apply(fun, [el]) for el in lis]], []))
    #uniquify(sum([[apply(fun,[el]) for el in lis]],[]))

def uniquify(lis):
    dic = {}
    for el in lis:
        dic[str(el)] = 1
    return [eval(k) for k in dic.keys()]

def parse_prefix(stem_code):
    '''
    Identifies whether a stem starts with 'm'+vowel or 't'+vowel 
    and augments the (prefix,stem,code) options accordingly.
    
    Parameters:
    ----------------
    stem_code is a tuple pair (stem, code) of strings
    
    Returns:
    -----------
    A list of triples (prefix,stem,code) of strings, with augmented codes.
    
    Examples:
    -------------
    >>> print BANARM.parse_prefix(('muriyH',''))
    >>> [('', 'muriyH', ''), ('mu', 'riyH', 'PFX[m] ')]
    '''
    (stem,code) = stem_code
    prefix_stem_code_list = [('',stem,code)]
    if len(stem)<2:
        return prefix_stem_code_list
    else:
        if stem[0] in 'tm' and stem[1] in vowels:
            prefix_stem_code_list.append((stem[0:2],stem[2:],code+'PFX['+stem[0]+'] '))
    return prefix_stem_code_list

def substitute_initial(prefix_stem_code_list):
    '''
    Examples:
    -------------
    >>>  print(BANARM.substitute_initial([('','>us~','')]))
    >>>  [('', '>us~', ''), ('', 'Aus~', 'IN '), ('', "'us~", 'IN '), ('', '<us~', 'IN '), ('', '{us~', 'IN ')]
    '''
    initialed_prefix_stem_code_list=[]
    for (prefix,stem,code) in prefix_stem_code_list:
        initialed_prefix_stem_code_list.append((prefix,stem,code))
        if len(stem)>0 and stem[0] in initials:
            for init in initials:
                if init != stem[0]:
                    initialed_prefix_stem_code_list.append((prefix,init+stem[1:],code+'IN '))
    return initialed_prefix_stem_code_list

def geminate_tildes(prefix_stem_code):
    '''
    Examples:
    -------------
    >>>  print(BANARM.geminate_tildes([('','>us~','')]))
    >>>  [('', '>us~', ''), ('', '>uss', 'GM/3 ')]
    '''
    result = [prefix_stem_code]
    prefix, stem, code = prefix_stem_code
    tildes = find_tilde_locations(stem)     
    for p in powerset_graycode(tildes):
        tmp = []
        pcode=''
        for i in range(len(stem)):
            if i not in p:
                tmp += [stem[i]]
            elif i in p and i>0:
                tmp += ['a', stem[i-1]]
                pcode=pcode+'/'+str(i+1)
        if len(pcode)>0:
            result += [(prefix, string.join(tmp,''), code+'GM'+pcode+' ')]
    return result

def find_tilde_locations(stem):
    '''
    List places where stem has a tilde
    
    Examples:
    -------------
    >>>  print(BANARM.find_tilde_locations('abc~de~'))
    >>>  [3, 6]

    '''
    return [index for (index,char) in enumerate(stem) if char == '~']

def powerset_graycode(s):
    '''
    Helper function.
    Creates a generator object for all subsets from a list or set s.
    INS can't remember where he found this code.
    Perhaps this should go in a combinatorics utility module?
    
    Example:
    ------------
    >>> print([list(places) for places in BANARM.powerset_graycode([3,6])])
    >>> [[], [3], [3, 6], [6]]
    '''
    d = dict(zip(
            (1<<i for i in range(len(s))),
            (set([e]) for e in s)
            ))
    subset = set()
    yield list(subset)
    for i in range(1, 1<<len(s)):
        subset = subset ^ d[i & -i]
        yield list(subset)

def two_consonant_expand(prefix_stem_code_list):
    '''
    Examples:
    -------------
    >>>  print(BANARM.two_consonant_expand([('','>us~','')]))
    >>>  [[('', '>us~', '')], [('>ius~', 'CC[i] '), ('>yus~', 'CC[y] '), ('>uss~', 'CC2 ')]]

    '''
    expanded_prefix_stem_code_list = []
    for (prefix,stem,code) in prefix_stem_code_list:
        expanded_prefix_stem_code_list.append([(prefix,stem,code)])
        places = find_root_locations(stem)
        if len(places) == 2 and places[1]-places[0] == 2:
            expanded_prefix_stem_code_list.append([
                (prefix,stem[:places[0]+1]+"i"+stem[places[1]-1:],code+'CC[i] '),
                (prefix,stem[:places[0]+1]+"y"+stem[places[1]-1:],code+'CC[y] '),
                (prefix,stem[:places[1]]+stem[places[1]]+stem[places[1]:],code+'CC2 ')
                ])
    return [inner for outer in expanded_prefix_stem_code_list for inner in outer]
  
def find_root_locations(stem):
    """List places where stem has potential root""" 
    return [x for x in range(len(stem)) if stem[x] not in 'aeiou~']
    
def find_nonvowel_locations(stem):
    """List places where stem has non-vowels""" 
    return [x for x in range(len(stem)) if stem[x] not in 'aeiou~Awy']

def meta_expand(prefix_stem_code):
    '''
    >>>  print(BANARM.meta_expand([('','DAfiy','')]))
    >>>  set([('', 'D%fiy', 'V/5 '), ('', 'D%fi#', 'C/5 ')])
    '''
    expanded = set([])
    prefix,stem,code = prefix_stem_code
    stem=add_meta(stem)
    places = find_meta_characters(stem)
    for g in powerset_graycode(places):
        gcode='C/'
        for gp in g:
            gcode=gcode+str(gp+1)
        if len(g) == 0:
            gcode = ''
        h = places.difference(g)
        hcode='V/'
        for hp in h:
            hcode=hcode+str(hp+1)
        if len(h) == 0:
            hcode = ''
        expanded.add((prefix, list_meta_expand(stem,g,h),
                      code+gcode+hcode+' '))
    return list(expanded)

def consonant_meta(stem):
    """Intermediate meta-characters [] are replaced by [@#=]"""
    stem = stem.replace('%','A')
    stem = stem.replace('#','y')
    stem = stem.replace('=','w')
    return stem

def meta_expand_list(psc_list):
    return sum([meta_expand(outer) for outer in psc_list],[])
        
def add_meta(stem):
    '''
    Helper function for meta_expand function
    [Ayw] are replaced by intermediate meta-characters [%+?]
    '''
    
    stem = stem.replace('A','%')
    stem = stem.replace('y','+')
    stem = stem.replace('w','?')
    return stem 

def find_meta_characters(stem):
    '''
    Helper function for meta_expand function
    Does what it say on the tin
    '''
    
    return set([x for x in range(len(stem)) if stem[x] in '&+?'])

def list_meta_expand(stem,placesx,placesy):
    '''
    Helper function for meta_expand function
    Does what it say on the tin
    '''
    
    s = list(stem)
    ''' 
        meta_characters -> consonants
    '''
    for p in placesx:
        if s[p] == '%':
            s[p] = '@'    # = -> A as consonant
        elif s[p] == '+':
            s[p] = '#'    # # -> y
        elif s[p] == '?':
            s[p] = '='    # = -> w
    '''
        meta_characters -> long vowels
    '''
    for p in placesy:
        if s[p] == '%':
            s[p] = 'A'
        elif s[p] == '+':
            s[p] = 'y'
        elif s[p] == '?':
            s[p] = 'w'
    lme = ''.join(s) 
    return lme

def make_parse_dictionary():
    '''
    Creates a python dictionary structure from the corpus in the file 
    annotated-dictionary.xls, omitting nonwords of various types. 
    The keys are valid the stems and the values are the corresponding 
    parsed (root, pattern) pairs in the Dictionary.
    '''
    nonwordtypes = ['abbreviation', 'compound', 'dialect_word', 'foreign_word', \
                'function_word', 'interjection', 'letter_name', 'proper_name']

    baam_parsings=baampath + 'annotated-dictionary.xls'
    parsings=form.loadexcel(baam_parsings)
    source_parsings=parsings.sheet_by_name('dictionary')
    entries=source_parsings.col_values(0)
    root=source_parsings.col_values(1)
    pattern=source_parsings.col_values(2)
    source_types=parsings.sheet_by_name('types')
    types = source_types.col_values(0)
    dictionarydic={}
    stemdic={}
    point=0
    for (point, entry) in enumerate(entries):
        if pattern[point] not in nonwordtypes:
            dictionarydic[entry]=(root[point],pattern[point])
    return dictionarydic


def make_BANARM_dictionary(source):
    ''' 
    Create a BANARM dictionary for a named tab ('source') in the
    spreadsheet 'root-pattern-frequencies.xls'.

    Requires module 'form'
    
    Parameter:
    ---------------
    source: string 
        Get BAAM data from Excel spreadsheet using function from module 'form'
    '''
    baam_databook=baampath+r'root-pattern-frequencies.xls'
    databook=form.loadexcel(baam_databook)
    source_sheet = databook.sheet_by_name(source)
    source_names = source_sheet.col_values(0)
    source_tokens = source_sheet.col_values(1)
    source_no = source_sheet.nrows
    return make_dictionary(source_names,source_tokens)

def make_dictionary(names,tokens):
    """Function to create a dictionary from lists of strings and frequencies"""
    dic = {}
    for (point, name) in enumerate(names):
        dic[name] = tokens[point]
    return dic

def dictionary_sets(stems,dictionary,report=True):
    '''Report statistics for relationship between corpus and dictionary'''
    dictionary_stems = set(dictionary.keys())
    corpus_stems = set(stems.keys())
    unique_dictionary_stems = list(dictionary_stems.difference(corpus_stems))
    unique_corpus_stems = list(corpus_stems.difference(dictionary_stems))
    common_stems = list(dictionary_stems.intersection(corpus_stems))
    commentary = ''
    commentary +=  "There are %d stems in the dictionary\n" % (len(dictionary_stems))
    commentary += "There are %d stems in the corpus\n" % (len(corpus_stems))
    commentary += "There are %d stems unique to the dictionary\n" % (len(unique_dictionary_stems))
    commentary += "There are %d stems unique to the corpus\n" % (len(unique_corpus_stems))
    commentary += "There are %d stems common to both corpus and dictionary\n"%(len(common_stems))
    #print commentary
    return {'commentary': commentary, 'stems': dictionary_stems,
            'corpus': corpus_stems, 'common': common_stems}

#print analyst(('mu', 'riyH', 'PFX[m] '))

##    expanded_stem_code_list = two_consonant_expand(stem,code)
##    parts = set([])
##    for (expanded_stem,code) in expanded_stem_code_list:
##        for (meta_expanded_stem,code) in meta_expand(expanded_stem,code):
##            parts.update(parse(meta_expanded_stem,code))
##    return parts

def count(dic,entry):
    if entry in dic:
        return int(dic[entry])
    else:
        return 0

def gather(prpc_list):
    gathered = {}
    for prpc in prpc_list:
        prefix, root, pattern, code = prpc
        parsing = (root, prefix+pattern)
        if parsing in gathered:
            gathered[parsing].append(code)
        else:
            gathered[parsing] = [code]
    return gathered

'''Create parsed dictionary'''
parsed_dictionary = make_parse_dictionary()
#print len(parsed_dictionary), 'parsed in dictionary'

'''Create frequency dictionaries for stems, roots and patterns'''
stem_dictionary = make_BANARM_dictionary('stem_pointed_token')
#print len(stem_dictionary), 'corpus stems'
root_dictionary = make_BANARM_dictionary('root_token')
#print len(root_dictionary), 'corpus roots'
pattern_dictionary = make_BANARM_dictionary('pattern_token')
#print len(pattern_dictionary), 'corpus patterns'
#print 'making dictionary sets ...'
dic_sets = dictionary_sets(stem_dictionary, parsed_dictionary, True)
print dic_sets['commentary']

test_cases =  [dic_sets['common'][i] for i in [1000,2000,3000]]

def evaluate(stem):
    options = analyst(('',stem,''))
    g_options = gather(options)
    solutions = []
    valid = 0
    invalid = 0
    for rp in g_options.keys():
        #print prsc
        codes =  g_options[rp]
        root, pattern = rp
        root_count = count(root_dictionary, root)
        pattern_count = count(pattern_dictionary, pattern)
        if root_count * pattern_count > 0:
            solutions += (stem, (root, pattern), root, root_count, pattern, pattern_count, codes)
            valid += 1
        else:
            invalid += 1
    if stem in parsed_dictionary:
        dic_sol = parsed_dictionary[stem]
    else:
        dic_sol = None
    return {'stem': stem, 'solutions': solutions, 'valid': valid, 'invalid': invalid, 'dictionary': dic_sol}

solutions={}
solfile = open('solutions','wb')
for stem in test_cases:
    solutions[stem]=evaluate(stem)
cPickle.dump(solutions,solfile)
print solutions



#print root, count(root_dictionary,root), pattern, count(pattern_dictionary,pattern)    

##
##def print_entry(file, stem, prefix_body_root_pattern_code_list):
##    for (prefix,body,root,rootcount,pattern,patterncount,comment,code) in prefix_body_root_pattern_code_list:
##        if prefix != '':
##            print >> file, '%s: [%s]%s ... %s (%s, %s) [%d, %d] %s' % (stem,prefix,body,comment,root,pattern,rootcount,patterncount,code)
##        else:
##            print >> file, '%s:   %s ... %s (%s, %s) [%d, %d] %s' % (stem,body,comment,root,pattern,rootcount,patterncount,code)
##
##def print_solutions(file, stem, solutions):
##            for (root,pattern,logfreq) in solutions:
##                print >> file, '%s: (%s, %s) [%d]' % (stem,root,pattern,logfreq)
##                        
##def parse(stem,code):
##    raw_code=code
##    parts = set([])
##    stemlist=list(stem)
##    places=find_nonvowel_locations(stem)
##    if len(places) >= 3:
##        if '3R' not in raw_code:
##            code=raw_code+'3R '
##        for root in generate_combinations(places,3):
##            try_pattern=list(consonant_meta(stem))
##            try_root=string.join([try_pattern[x] for x in root],'')
##            for x in range(3):
##                try_pattern[root[x]]=slot_names3[x]
##            try_pattern=string.join(try_pattern,'')
####            if try_pattern[-1] in 'oN':
####                parts.add((try_root,try_pattern[:-1],code+'PattTrunc[oN.]'))
##            if try_root[-1] in 'wyY':
##                for char in 'wyY':
##                    if try_root[-1] != char:
##                        parts.add((try_root[:-1]+char,try_pattern,code+'RootEnd[wyY] '))
##            parts.add((try_root,try_pattern,code))
##            if try_pattern[-1] in 'wyY':
##                for char in 'wyY':
##                    if try_pattern[-1] != char:
##                        parts.add((try_root,try_pattern[:-1]+char,code+'PattEnd[wyY] '))
##            if try_pattern[-1] in "&'}":
##                for char in "&'}":
##                    parts.add((try_root,try_pattern[:-1]+char,code+'PattEnd[&] '))
####            if try_pattern[:2] in ["{i", "<i"]:
####                for first_two  in ["{i", "<i"]:
####                    if first_two != try_pattern[:2]:
####                        parts.add((try_root,first_two+try_pattern[2:],code+'PattInit[{i |<i] '))
##            if try_pattern[0] in initials:
##                for first in initials:
##                    if first != try_pattern[0]:
##                        parts.add((try_root,first+try_pattern[1:],code+'PattInitGlot '))
##            if try_root[1] in glottal_medial_root:
##                for char in glottal_medial_root:
##                    if try_root[1] != char:
##                        parts.add((try_root[0]+char+try_root[2],try_pattern,code+'RootMedGlot '))             
##            if try_root[2] in glottal_final_root:
##                for char in glottal_final_root:
##                    if try_root[1] != char:
##                        parts.add((try_root[:2]+char,try_pattern,code+'RootFinGlot '))             
##    if len(places) >= 4:
##        if '4R' not in raw_code:
##            code=raw_code+'4R '
##        for root in generate_combinations(places,4):
##            try_pattern=list(consonant_meta(stem))
##            try_root=string.join([try_pattern[x] for x in root],'')
##            for x in range(4):
##                try_pattern[root[x]]=slot_names4[x]
##            try_pattern=string.join(try_pattern,'')
####            if try_pattern[-1] in 'oN':
####                parts.add((try_root,try_pattern[:-1],code+'FP-'+try_pattern[-1]+' '))
####            else:
##            parts.add((try_root,try_pattern,code))
##            if try_pattern[-1] in "yY":
##                for char in "yY":
##                    if try_pattern[1] != char:
##                        parts.add((try_root,try_pattern[:-1]+char,code+'PattEnd[y|Y] '))
##            if try_pattern[-1] in "&'}":
##                for char in "&'}":
##                    if try_pattern[1] != char:
##                        parts.add((try_root,try_pattern[:-1]+char,code+"PattEnd[&'}] "))
####            if try_pattern[:2] in ['{i', '<i']:
####                for first_two  in ['{i', '<i']:
####                    if first_two != try_pattern[:2]:
####                        parts.add((try_root,first_two+try_pattern[2:],code+'PattEnd[{i |<i] '))
##            if try_pattern[0] in initials:
##                for first in initials:
##                    if first != try_pattern[0]:
##                        parts.add((try_root,first+try_pattern[1:],code+'PattInitGlot '))
##    cvc = set([])
##    for (try_root,try_pattern,code) in parts:
##        last=try_pattern.find('E')
##        if try_root[1] == 'A':
##            cvc.add((try_root[0]+"y"+try_root[2],try_pattern[:last]+"aEa"+
##                try_pattern[last+1:],code+'Sub[A/(y,aEa)] '))
##            cvc.add((try_root[0]+"w"+try_root[2],try_pattern[:last]+"aEa"+
##                try_pattern[last+1:],code+'Sub[A/(w,aEa)] '))
##        if try_root[1] == 'w':
##            cvc.add((try_root[0]+"w"+try_root[2],try_pattern[:last]+"uEu"+
##                try_pattern[last+1:],code+'Sub[w/(w,uEu)] '))
##        if try_root[1] == 'y':
##            cvc.add((try_root[0]+"y"+try_root[2],try_pattern[:last]+"iEi"+
##                try_pattern[last+1:],code+'Sub[y/(y,iEi)] '))
##    parts.update(cvc)
####    oddsandends = set([])
####    for (root,pattern) in parts:
####        if root[2] in "Yyw":
####            for ending in "Yyw":
####                oddsandends.add((root[0:2]+ending,pattern))
####        if root[1] in "Awy":
####            for middle in "Awy":
####                oddsandends.add((root[0]+middle+root[2],pattern))
####        if root[2] in glottal_stops:
####            for ending in glottal_stops:
####                oddsandends.add((root[0:2]+ending,pattern))
####    parts.update(oddsandends)
##    return parts
##
##def find_nonvowel_locations(stem):
##    """List places where stem has non-vowels""" 
##    return [x for x in range(len(stem)) if stem[x] not in 'aeiou~Awy']
##
##def generate_combinations(items, n):
##    """Combine sets of size n from items
##    
##    Perhaps this should go in a combinatorics utility module?"""
##    if n == 0:
##        yield []
##    else:
##        for i in xrange(len(items)):
##            for cc in generate_combinations(items[i+1:], n-1):
##                yield [items[i]] + cc
##
##def consonant_meta(stem):
##    """Intermediate meta-characters [] are replaced by [@#=]"""
##    stem = stem.replace('%','A')
##    stem = stem.replace('#','y')
##    stem = stem.replace('=','w')
##    return stem
##
##def root_pattern_match(r0,p0,r1,p1,strictness='strict'):
##    if r0 == r1:
##        if p0 == p1:
##            return True
##        elif strictness == 'strict':
##            return p0[-1] == 'a' and p0[:-1] ==p1
##        elif strictness == 'EaEoEua':
##            re_string1 = re.compile('(Ea|Eo|Eu)')
##            re_string2 = re.compile('a')
##            return re_string2.sub('',re_string1.sub('E',p0)) == re_string2.sub('',re_string1.sub('E',p1))
##        else:
##            return False
##    else:
##        return False
##
##'''Report statistics for relationship between corpus and dictionary'''
##dictionary_stems, corpus_stems, common_stems = dictionary_sets(stem_dictionary, parsed_dictionary)
##
##def make_full_analyses(dictionary=stem_dictionary, limit=100000):
##
##    analyses = {}
##
##    #for (count, raw_stem) in enumerate(common_stems):
##    for (count, raw_stem) in enumerate(dictionary):
##        if count > limit:
##            break
##        else:
##            #print count, raw_stem
##            analyses[raw_stem]={'parsings':[], 'solutions':set([])}
##            for (prefix,body) in geminate_tildes(substitute_initial(parse_prefix(raw_stem))):
##                for (root,pattern) in analyse(body):
##                    if root_dictionary.has_key(root):
##                        rootcount = root_dictionary[root]
##                    else:
##                        rootcount = 0
##                    if pattern_dictionary.has_key(prefix+pattern):
##                        patterncount = pattern_dictionary[prefix+pattern]
##                    else:
##                        patterncount = 0
##                    if rootcount*patterncount==0:
##                        comment = ''
##                    else:
##                        comment = '*'
##                        analyses[raw_stem]['solutions'].add((root,prefix+pattern,10*np.log2(rootcount*patterncount)))
##                    analyses[raw_stem]['parsings'].append((prefix,body,root,rootcount,prefix+pattern,patterncount,comment))
##                    
##        #print_entry(logfile, raw_stem, analyses[raw_stem]['parsings'])
##    return analyses
##
##def report_detailed_analyses(dictionary=stem_dictionary, limit=100000):
##
##    analyses = {}
##
##    #for (count, raw_stem) in enumerate(common_stems):
##    for (count, raw_stem) in enumerate(dictionary):
##        if count > limit:
##            break
##        else:
##            #print count, raw_stem
##            analyses[raw_stem]={'parsings':[], 'solutions':set([])}
##            for (prefix,body,code) in geminate_tildes(substitute_initial(parse_prefix(raw_stem,''))):
##                for (root,pattern,code) in analyse(body,code):
##                    if root_dictionary.has_key(root):
##                        rootcount = root_dictionary[root]
##                    else:
##                        rootcount = 0
##                    if pattern_dictionary.has_key(prefix+pattern):
##                        patterncount = pattern_dictionary[prefix+pattern]
##                    else:
##                        patterncount = 0
##                    if rootcount*patterncount==0:
##                        comment = ' '
##                    else:
##                        comment = '*'
##                        analyses[raw_stem]['solutions'].add((root,prefix+pattern,10*np.log2(rootcount*patterncount),code))
##                    analyses[raw_stem]['parsings'].append((prefix,body,root,rootcount,prefix+pattern,patterncount,comment,code))
##                    
##        print_entry(logfile, raw_stem, analyses[raw_stem]['parsings'])
##        
##        dictionarydic = make_parse_dictionary()
##        if raw_stem in common_stems:
##            (dic_root, dic_pattern) = dictionarydic[raw_stem]
##            dic_root = str(dic_root)
##            dic_pattern = str(dic_pattern)
##            print >> logfile,'Dic      '+raw_stem+':     ('+dic_root+', '+dic_pattern+')\n'
##        else:
##            print >> logfile,'No dictionary entry for ',raw_stem,'\n'
##
##    return analyses
##
##def make_terse_analyses(dictionary=stem_dictionary, limit=100000):
##
##    analyses = {}
##
##    #for (count, raw_stem) in enumerate(common_stems):
##    for (count, raw_stem) in enumerate(dictionary):
##        if count > limit:
##            break
##        else:
##            #print count, raw_stem
##            analyses[raw_stem]={'parsings':[], 'solutions':set([])}
##            for (prefix,body) in geminate_tildes(substitute_initial(parse_prefix(raw_stem))):
##                for (root,pattern) in analyse(body):
##                    if root_dictionary.has_key(root):
##                        rootcount = root_dictionary[root]
##                    else:
##                        rootcount = 0
##                    if pattern_dictionary.has_key(prefix+pattern):
##                        patterncount = pattern_dictionary[prefix+pattern]
##                    else:
##                        patterncount = 0
##                    if rootcount*patterncount==0:
##                        comment = ''
##                    else:
##                        comment = '*'
##                        analyses[raw_stem]['solutions'].add((root,prefix+pattern,10*np.log2(rootcount*patterncount)))
##                    #analyses[raw_stem]['parsings'].append((prefix,body,root,rootcount,prefix+pattern,patterncount,comment))
##                    
##        #print_entry(logfile, raw_stem, analyses[raw_stem]['parsings'])
##    return analyses
##
##logfile = open(baampath+'BAAM-09-02-10-solutions.log', 'wt')
##
###analyses = make_full_analyses(dictionary=common_stems, limit=len(common_stems))
###testlist = ['muriyH','>us~','sum~','Tuwl','$uwraY']
##testlist = common_stems
##analyses = report_detailed_analyses(dictionary=testlist, limit = 250 )
##
###analyses = report_detailed_analyses(dictionary=['muriyH'])
##
####solutions_by_counts = {}
####
####for key in analyses.keys():
####    try:
####        number_of_solutions = len(analyses[key]['solutions'])
####    except:
####        number_of_solutions = 0
####    try:
####        solutions_by_counts[number_of_solutions].append(key)
####    except:
####        solutions_by_counts[number_of_solutions] = [key]
####
####for count in solutions_by_counts.keys():
####    print >> logfile, count, len(solutions_by_counts[count])
####
####dictionarydic = make_parse_dictionary()
####
####re_string1 = re.compile('(Ea|Eo|Eu)')
####re_string2 = re.compile('(a)')
####    
####total_listed = 0
####total_unlisted = 0
####total_agree = 0
####total_disagree = 0
####
####for count in solutions_by_counts.keys():
####
####    listed = 0
####    unlisted = 0
####    agree = 0
####    disagree = 0
####
####    if count == 0:
####        print >> logfile, 'No solutions found ...'
####        for unsolved in solutions_by_counts[0]:
####            if unsolved in common_stems:
####                (dic_root, dic_pattern) = dictionarydic[unsolved]
####                dic_root = str(dic_root)
####                dic_pattern = str(dic_pattern)
####                print >> logfile,unsolved ,': (', dic_root,', ',dic_pattern,')'
####                #print_entry(logfile, unsolved, analyses[unsolved]['parsings'])
####                listed += 1
####            else:
####                unlisted += 1
####    elif count == 1:
####        print >> logfile, 'One solution found ...'
####        for singleton in solutions_by_counts[1]:
####            (root,pattern,logfreq) = list(analyses[singleton]['solutions'])[0]
####            if singleton in common_stems:
####                (dic_root, dic_pattern) = dictionarydic[singleton]
####                dic_root = str(dic_root)
####                dic_pattern = str(dic_pattern)
####                #print (root,pattern,dic_root, dic_pattern,root==dic_root,re_string1.sub('E',pattern)==re_string1.sub('E',dic_pattern),root_pattern_match(root,pattern,dic_root, dic_pattern,strictness='EaEo'))
####                if root_pattern_match(root,pattern,dic_root, dic_pattern,strictness='EaEoEua'):
####                    agree += 1
####                else:
####                    print >> logfile, '[', singleton, '~ (', root,', ', pattern,') # (', dic_root,', ', dic_pattern,')]'
####                    #print_entry(logfile, singleton, analyses[singleton]['parsings'])
####                    disagree += 1
####    else:
####        print >> logfile, 'Several solutions found ... [',count,']'
####        for multiple_solutions in solutions_by_counts[count]:
####            if multiple_solutions in common_stems:
####                solution_list = list(analyses[multiple_solutions]['solutions'])
####                where = np.argmax([l for (r,p,l) in solution_list])
####                (root,pattern,logfreq) = solution_list[where]
####                (dic_root, dic_pattern) = dictionarydic[multiple_solutions]
####                dic_root = str(dic_root)
####                dic_pattern = str(dic_pattern)
####                if root_pattern_match(root,pattern,dic_root, dic_pattern,strictness='EaEoEu'):
####                    agree += 1
####                else:
####                    print >> logfile, '[', multiple_solutions, '~ (', root,', ', pattern,') # (', dic_root,', ', dic_pattern,')]'
####                    #print_solutions(logfile, multiple_solutions, analyses[multiple_solutions]['solutions'])
####                    disagree += 1
####    if count>0:
####        print (count,len(solutions_by_counts[count]), agree, disagree,(agree*100.)/(agree+disagree))
####    if count==0:
####        print (count,len(solutions_by_counts[count]), listed, unlisted)
####
####    total_listed += listed
####    total_unlisted += unlisted
####    total_agree  +=  agree
####    total_disagree += disagree
####
####
####print >> logfile, total_listed, 'Listed, ', total_unlisted, ' Unlisted'
####print >> logfile, total_agree, 'Agree, ', total_disagree, ' Disagree'
####print >> logfile, (total_agree*100.)/(total_agree+total_disagree), '% agreement'
#### 
#####anafile = open('analyses.pkl','wb')
#####cPickle.dump(analyses,anafile)
####
#####solfile = open('solutions_by_counts.pkl','wb')
#####cPickle.dump(solutions_by_counts,solfile)
####
#####anafile.close()
#####solfile.close()
##        
##logfile.close()

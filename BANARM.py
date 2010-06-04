'''
Author: Ian Nimmo-Smith
Description: Python library for Bayesian ANalysis of ARabic Morphology (BANARM)
Date: 04 May 2010
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

def unique_strings(lis):
    dic = {}
    for el in lis:
        dic[el] = 1
    return [k for k in dic.keys()]

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
            dictionarydic[entry]=(root[point].encode(),pattern[point].encode())
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
    for (n, name) in enumerate(source_names):
        source_names[n] = name.encode()
    source_tokens = source_sheet.col_values(1)
    source_no = source_sheet.nrows
    return make_dictionary(source_names,source_tokens)

def make_dictionary(names,tokens):
    """Function to create a dictionary from lists of strings and frequencies"""
    dic = {}
    for (point, name) in enumerate(names):
        dic[name] = tokens[point]
    return dic

def dictionary_sets(stems,dictionary,report=True,save=False):
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
    if report:
        print commentary
    stem_sets = {'commentary': commentary, 'stems': dictionary_stems,
            'corpus': corpus_stems, 'common': common_stems}
    if save:
        stem_sets_file = open('stem_sets','wb')
        cPickle.dump(stem_sets,stem_sets_file)
    return stem_sets

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

def evaluate(stem,root_dictionary, pattern_dictionary, parsed_dictionary):
    stem=stem.encode()
    options = analyst(('',stem,''))
    g_options = gather(options)
    solutions = []
    roots = []
    patterns = []
    valid = 0
    invalid = 0
    for rp in g_options.keys():
        #print prsc
        codes =  g_options[rp]
        root, pattern = rp
        root_count = count(root_dictionary, root)
        if root_count > 0:
            roots.append((root, root_count))
        pattern_count = count(pattern_dictionary, pattern)
        if pattern_count > 0:
            patterns.append((pattern, pattern_count))
        if root_count * pattern_count > 0:
            solutions += (stem, (root, pattern), root, root_count, pattern, pattern_count, codes)
            valid += 1
        else:
            invalid += 1
    #print roots
    roots = unique_strings(roots)
    #print patterns
    patterns = unique_strings(patterns)
    if stem in parsed_dictionary:
        dic_sol = parsed_dictionary[stem]
    else:
        dic_sol = None
    return {'stem': stem, 'roots': roots, 'patterns': patterns,'solutions': solutions, 'valid': valid, 'invalid': invalid, 'dictionary': dic_sol}


def create_solutions():
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

    #test_cases =  [dic_sets['common'][i] for i in [1000,2000,3000]]

    solutions={}
    solfile = open('fullsolutions','wb')
    solcount = 0

    for stem in stem_dictionary:
        solcount += 1
        if solcount % 1000 == 0:
            print '... entry ',solcount
        solutions[stem.encode()] = evaluate(stem, root_dictionary, pattern_dictionary, parsed_dictionary)

    # this took 6 hours on devel06 on 30-Apr-2010
    # next version took 6.5 h on 2 May 2010

    cPickle.dump(solutions,solfile)
    solfile.close()
    return solutions

def load_solutions():
    print '... starting loading solutions ...'
    solfile = open('fullsolutions','r')
    solutions= cPickle.load(solfile)
    print '... finishing loading ...'
    return solutions

def load_stem_sets():
    print '... starting loading stem sets ...'
    stem_sets_file = open('stem_sets','r')
    stem_sets= cPickle.load(stem_sets_file)
    print '... finishing loading ...'
    return stem_sets

def root_pattern_match(rp_a,rp_b):
    r_a, p_a = rp_a
    r_b, p_b = rp_b
    if p_a[-1] in 'Na':
        p_a = p_a[:-1]
    if p_b[-1] in 'Na':
        p_b = p_b[:-1]
    if p_a == p_b:
        return True
    else:
        return False

def maximum_frequency_root(solution):
    r_f = solution['roots']
    if len(r_f) == 0:
        return ('',0)
    else:
        return r_f[np.argmax([freq for (root,freq) in r_f])]
    
def maximum_frequency_pattern(solution):
    p_f = solution['patterns']
    if len(p_f) == 0:
        return ('',0)
    else:
        return p_f[np.argmax([freq for (pattern,freq) in p_f])]
    
if __name__ == '__main__':

    make_dics = True

    if make_dics:
        '''Create parsed dictionary'''
        parsed_dictionary = make_parse_dictionary()
        '''Create frequency dictionaries for stems, roots and patterns'''
        stem_dictionary = make_BANARM_dictionary('stem_pointed_token')
        root_dictionary = make_BANARM_dictionary('root_token')
        pattern_dictionary = make_BANARM_dictionary('pattern_token')
        stem_sets = dictionary_sets(stem_dictionary, parsed_dictionary, report=True, save=True)
    else:
        stem_sets =  load_stem_sets()

    common = stem_sets['common']

    import time
    starttime = time.time()

    #solutions = create_solutions()
    solutions = load_solutions()

    none = 0
    single = 0
    multiple = 0
    solvedcount = 0
    unsolvedcount = 0

    unsolved = open('unsolved','w')
    invalid = open('invalid', 'w')

    maxroot_hit = 0
    maxroot_miss = 0

    maxpattern_hit = 0
    maxpattern_miss = 0

    multmaxroot_hit = 0
    multmaxroot_miss = 0

    multmaxpattern_hit = 0
    multmaxpattern_miss = 0

    #for stem in [common[j] for j in range(401,10000,500)]:
    for stem in common:
        dictionary_solution = parsed_dictionary[stem]
        #print dictionary_solution
        s = solutions[stem]
        #print maximum_frequency_root(s)
        if dictionary_solution[0] == maximum_frequency_root(s)[0]:
            maxroot_hit +=1
        else:
            maxroot_miss +=1
        if dictionary_solution[1] == maximum_frequency_pattern(s)[0]:
            maxpattern_hit +=1
        else:
            maxpattern_miss +=1
        if s['valid'] == 0:
            print >> invalid, s['stem'], s['dictionary']
            none += 1
        elif s['valid'] == 1:
            single +=1
            if root_pattern_match(s['solutions'][1], s['dictionary']):
                solvedcount += 1
            else:
                unsolvedcount += 1
                print >> unsolved, s['stem'], s['solutions'][1], '!=', s['dictionary']
        else:
            multiple += 1
            if dictionary_solution[0] == maximum_frequency_root(s)[0]:
                multmaxroot_hit +=1
            else:
                multmaxroot_miss +=1
            if dictionary_solution[1] == maximum_frequency_pattern(s)[0]:
                multmaxpattern_hit +=1
            else:
                multmaxpattern_miss +=1
            maxpattern_miss +=1
    unsolved.close()
    invalid.close()
    print 'none %d\n \
           single %d\n \
           solved %d\n \
           unsolved %d\n \
           multiple solutions %d\n \
           maxroot_hits %d\n \
           maxroot_miss %d\n \
           maxpattern_hits %d\n \
           maxpattern_miss %d\n \
           multnmaxroot_hits %d\n \
           multmaxroot_miss %d\n \
           multmaxpattern_hits %d\n \
           multmaxpattern_miss %d' \
           % (none,single,solvedcount,unsolvedcount,multiple,maxroot_hit, \
               maxroot_miss, maxpattern_hit, maxpattern_miss, multmaxroot_hit, multmaxroot_miss, multmaxpattern_hit, multmaxpattern_miss) 
    endtime = time.time()
    print '... program took %d seconds ' % np.int(endtime-starttime)



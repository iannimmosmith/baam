'''
Author: Ian Nimmo-Smith
Description: Python library for Bayesian Analysis of Arabic Morphology (BAAM)
Date: 09 Feb 2010
'''
import form
import string
import cPickle
import numpy as np
import re

slot_names3=['f', 'E', 'l']
slot_names4=['f', 'E', 'l', 'l']
initials = "A'<>{"
glottal_stops = "'|>&<}`{"
glottal_final_root = "'>A<{}"
glottal_medial_root = "&A'>"
extra_letters = "'Atslmnwyo"
pointing = 'aiu~oN'
vowels = 'aeiou~Awy'

baampath = r'/home/ian/Sami/baam/'

def make_parse_dictionary():
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
            if pattern[point][-1] in 'aN':
                dictionarydic[entry]=(root[point],pattern[point][:-1])
            else:
                dictionarydic[entry]=(root[point],pattern[point])
    return dictionarydic

def make_baam_dictionary(source):
    ''' Create a BAAM dictionary for a named tab ('source')
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
    dictionary_stems = keyset(dictionary)
    corpus_stems = keyset(stems)
    unique_dictionary_stems = dictionary_stems.difference(corpus_stems)
    unique_corpus_stems = corpus_stems.difference(dictionary_stems)
    common_stems = dictionary_stems.intersection(corpus_stems)
    if report:
        commentary = open(baampath+'BAAM.commentary', 'wt')
        print >> commentary, "There are %d stems in the dictionary" % len(dictionary_stems)
        print >> commentary, "There are %d stems in the corpus" % len(corpus_stems)
        print >> commentary, "There are %d stems unique to the dictionary" % len(unique_dictionary_stems)   
        print >> commentary, "There are %d stems unique to the corpus" % len(unique_corpus_stems)   
        print >> commentary, "There are %d stems common to both corpus and dictionary" % len(common_stems)
        commentary.close()
    return dictionary_stems, corpus_stems, common_stems
    
def keyset(dic):
    return set(dic.keys())

def parse_prefix(raw_stem,code):
    prefix_body_code_list = [('',raw_stem,code)]
    if len(raw_stem)<2:
        return prefix_body_list
    else:
        if raw_stem[0] in 'tm' and raw_stem[1] in vowels:
            prefix_body_code_list.append((raw_stem[0:2],raw_stem[2:],code+'Prefix['+raw_stem[0]+'] '))
    return prefix_body_code_list

def analyse(stem,code):
##    expanded_stem_list = []
##    interim_list = two_consonant_expand(stem)
##    for interim_stem in interim_list:
##        expanded_stem_list.extend(glottal_swap(stem))
##    parts = set([])
##    for expanded_stem in expanded_stem_list:
##        for meta_expanded_stem in meta_expand(expanded_stem):
##            parts.update(parse(meta_expanded_stem))
##    return parts
    expanded_stem_code_list = two_consonant_expand(stem,code)
    parts = set([])
    for (expanded_stem,code) in expanded_stem_code_list:
        for (meta_expanded_stem,code) in meta_expand(expanded_stem,code):
            parts.update(parse(meta_expanded_stem,code))
    return parts

def two_consonant_expand(stem,code):
    places = find_root_locations(stem)
    if len(places) == 2 and places[1]-places[0] == 2:
        expanded = [(stem[:places[0]+1]+"i"+stem[places[1]-1:],code+'2Cons[i] '),
        #expanded = [(stem[:places[0]+1]+"w"+stem[places[1]-1:],code+'CCw.'),
                    (stem[:places[0]+1]+"y"+stem[places[1]-1:],code+'2Cons[y] '),
                    (stem[:places[1]]+stem[places[1]]+stem[places[1]:],code+'2Cons2 ')]
        return expanded
    else:
        return [(stem,code)]

def find_root_locations(stem):
    """List places where stem has potential root""" 
    return [x for x in range(len(stem)) if stem[x] not in 'aeiou~']
    
'''Create parsed dictionary'''
parsed_dictionary = make_parse_dictionary()

'''Create frequency dictionaries for stems, roots and patterns'''
stem_dictionary = make_baam_dictionary('stem_pointed_token')
root_dictionary = make_baam_dictionary('root_token')
pattern_dictionary = make_baam_dictionary('pattern_token')

def print_entry(file, stem, prefix_body_root_pattern_code_list):
    for (prefix,body,root,rootcount,pattern,patterncount,comment,code) in prefix_body_root_pattern_code_list:
        if prefix != '':
            print >> file, '%s: [%s]%s ... %s (%s, %s) [%d, %d] %s' % (stem,prefix,body,comment,root,pattern,rootcount,patterncount,code)
        else:
            print >> file, '%s:   %s ... %s (%s, %s) [%d, %d] %s' % (stem,body,comment,root,pattern,rootcount,patterncount,code)

def print_solutions(file, stem, solutions):
            for (root,pattern,logfreq) in solutions:
                print >> file, '%s: (%s, %s) [%d]' % (stem,root,pattern,logfreq)
                        
def parse(stem,code):
    raw_code=code
    parts = set([])
    stemlist=list(stem)
    places=find_nonvowel_locations(stem)
    if len(places) >= 3:
        if '3R' not in raw_code:
            code=raw_code+'3R '
        for root in generate_combinations(places,3):
            try_pattern=list(consonant_meta(stem))
            try_root=string.join([try_pattern[x] for x in root],'')
            for x in range(3):
                try_pattern[root[x]]=slot_names3[x]
            try_pattern=string.join(try_pattern,'')
##            if try_pattern[-1] in 'oN':
##                parts.add((try_root,try_pattern[:-1],code+'PattTrunc[oN.]'))
            if try_root[-1] in 'wyY':
                for char in 'wyY':
                    if try_root[-1] != char:
                        parts.add((try_root[:-1]+char,try_pattern,code+'RootEnd[wyY] '))
            parts.add((try_root,try_pattern,code))
            if try_pattern[-1] in 'wyY':
                for char in 'wyY':
                    if try_pattern[-1] != char:
                        parts.add((try_root,try_pattern[:-1]+char,code+'PattEnd[wyY] '))
            if try_pattern[-1] in "&'}":
                for char in "&'}":
                    parts.add((try_root,try_pattern[:-1]+char,code+'PattEnd[&] '))
##            if try_pattern[:2] in ["{i", "<i"]:
##                for first_two  in ["{i", "<i"]:
##                    if first_two != try_pattern[:2]:
##                        parts.add((try_root,first_two+try_pattern[2:],code+'PattInit[{i |<i] '))
            if try_pattern[0] in initials:
                for first in initials:
                    if first != try_pattern[0]:
                        parts.add((try_root,first+try_pattern[1:],code+'PattInitGlot '))
            if try_root[1] in glottal_medial_root:
                for char in glottal_medial_root:
                    if try_root[1] != char:
                        parts.add((try_root[0]+char+try_root[2],try_pattern,code+'RootMedGlot '))             
            if try_root[2] in glottal_final_root:
                for char in glottal_final_root:
                    if try_root[1] != char:
                        parts.add((try_root[:2]+char,try_pattern,code+'RootFinGlot '))             
    if len(places) >= 4:
        if '4R' not in raw_code:
            code=raw_code+'4R '
        for root in generate_combinations(places,4):
            try_pattern=list(consonant_meta(stem))
            try_root=string.join([try_pattern[x] for x in root],'')
            for x in range(4):
                try_pattern[root[x]]=slot_names4[x]
            try_pattern=string.join(try_pattern,'')
##            if try_pattern[-1] in 'oN':
##                parts.add((try_root,try_pattern[:-1],code+'FP-'+try_pattern[-1]+' '))
##            else:
            parts.add((try_root,try_pattern,code))
            if try_pattern[-1] in "yY":
                for char in "yY":
                    if try_pattern[1] != char:
                        parts.add((try_root,try_pattern[:-1]+char,code+'PattEnd[y|Y] '))
            if try_pattern[-1] in "&'}":
                for char in "&'}":
                    if try_pattern[1] != char:
                        parts.add((try_root,try_pattern[:-1]+char,code+"PattEnd[&'}] "))
##            if try_pattern[:2] in ['{i', '<i']:
##                for first_two  in ['{i', '<i']:
##                    if first_two != try_pattern[:2]:
##                        parts.add((try_root,first_two+try_pattern[2:],code+'PattEnd[{i |<i] '))
            if try_pattern[0] in initials:
                for first in initials:
                    if first != try_pattern[0]:
                        parts.add((try_root,first+try_pattern[1:],code+'PattInitGlot '))
    cvc = set([])
    for (try_root,try_pattern,code) in parts:
        last=try_pattern.find('E')
        if try_root[1] == 'A':
            cvc.add((try_root[0]+"y"+try_root[2],try_pattern[:last]+"aEa"+
                try_pattern[last+1:],code+'Sub[A/(y,aEa)] '))
            cvc.add((try_root[0]+"w"+try_root[2],try_pattern[:last]+"aEa"+
                try_pattern[last+1:],code+'Sub[A/(w,aEa)] '))
        if try_root[1] == 'w':
            cvc.add((try_root[0]+"w"+try_root[2],try_pattern[:last]+"uEu"+
                try_pattern[last+1:],code+'Sub[w/(w,uEu)] '))
        if try_root[1] == 'y':
            cvc.add((try_root[0]+"y"+try_root[2],try_pattern[:last]+"iEi"+
                try_pattern[last+1:],code+'Sub[y/(y,iEi)] '))
    parts.update(cvc)
##    oddsandends = set([])
##    for (root,pattern) in parts:
##        if root[2] in "Yyw":
##            for ending in "Yyw":
##                oddsandends.add((root[0:2]+ending,pattern))
##        if root[1] in "Awy":
##            for middle in "Awy":
##                oddsandends.add((root[0]+middle+root[2],pattern))
##        if root[2] in glottal_stops:
##            for ending in glottal_stops:
##                oddsandends.add((root[0:2]+ending,pattern))
##    parts.update(oddsandends)
    return parts

def meta_expand(stem,code):
    expanded = set([])
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
        expanded.add((list_meta_expand(stem,g,h),code+gcode+hcode+' '))
    return expanded

def add_meta(stem):
    """[Ayw] are replaced by intermediate meta-characters [%+?]"""
    stem = stem.replace('A','%')
    stem = stem.replace('y','+')
    stem = stem.replace('w','?')
    return stem 

def find_meta_characters(stem):
    return set([x for x in range(len(stem)) if stem[x] in '&+?'])

def powerset_graycode(s):
    """Creates a generator object for all subsets from a set of size s
    
    Perhaps this should go in a combinatorics utility module?"""
    d = dict(zip(
            (1<<i for i in range(len(s))),
            (set([e]) for e in s)
            ))
    subset = set()
    yield subset
    for i in range(1, 1<<len(s)):
        subset = subset ^ d[i & -i]
        yield subset

def list_meta_expand(stem,placesx,placesy):
    s = list(stem)
    ''' 
        meta -> consonants
    '''
    for p in placesx:
        if s[p] == '%':
            s[p] = '@'    # = -> A as consonant
        elif s[p] == '+':
            s[p] = '#'    # # -> y
        elif s[p] == '?':
            s[p] = '='    # = -> w
    '''
        meta -> long vowels
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

def find_nonvowel_locations(stem):
    """List places where stem has non-vowels""" 
    return [x for x in range(len(stem)) if stem[x] not in 'aeiou~Awy']

def generate_combinations(items, n):
    """Combine sets of size n from items
    
    Perhaps this should go in a combinatorics utility module?"""
    if n == 0:
        yield []
    else:
        for i in xrange(len(items)):
            for cc in generate_combinations(items[i+1:], n-1):
                yield [items[i]] + cc

def consonant_meta(stem):
    """Intermediate meta-characters [] are replaced by [@#=]"""
    stem = stem.replace('%','A')
    stem = stem.replace('#','y')
    stem = stem.replace('=','w')
    return stem

def geminate_tildes(prefix_body_code_list):
    geminated_list=[]
##    tmp_list=[]
    for (prefix,body,code) in prefix_body_code_list:
        tildes = find_tilde_locations(body)     
        for p in powerset_graycode(tildes):
            tmp = list(body[:])
            pcode=''
            for i in p:
                if i>0:
                    tmp[i] = tmp[i-1]
                    pcode=pcode+'/'+str(i)
            if len(pcode)>0:
                geminated_list.append((prefix,string.join(tmp,''),code+'GM'+pcode+' '))
            else:
                geminated_list.append((prefix,string.join(tmp,''),code))
##            tmp_list.append((prefix,string.join(tmp,'')))
##    geminated_list=tmp_list[:]
##    for (prefix,body) in tmp_list:
##        if "t~" in body:
##            geminated_list.append((prefix,body.replace("t~","tw",1)))
    return geminated_list

def find_tilde_locations(stem):
    """List places where stem has a tilde""" 
    return [index for (index,char) in enumerate(stem) if char == '~']

##def glottal_swap(stem):
##    swapped = [stem]
##    for (i,gs) in enumerate(glottal_stops[:-1]):
##        if "A"+gs+"i" in stem:
##            for rgs in glottal_stops[i+1:]:
##                swapped.append(stem.replace("A"+gs+"i","A"+rgs+"i",1))
##    return swapped

def substitute_initial(prefix_body_code_list):              
    initialed_list=[]
    for (prefix,body,code) in prefix_body_code_list:
        if len(body)>0:
            if body[0] in initials:
                for init in initials:
                    if init != body[0]:
                        initialed_list.append((prefix,init+body[1:],code+'IN '))
                    else:
                        initialed_list.append((prefix,init+body[1:],code))
            else:
                initialed_list.append((prefix,body,code))
        else:
            initialed_list.append((prefix,body,code))
    return initialed_list

def root_pattern_match(r0,p0,r1,p1,strictness='strict'):
    if r0 == r1:
        if p0 == p1:
            return True
        elif strictness == 'strict':
            return p0[-1] == 'a' and p0[:-1] ==p1
        elif strictness == 'EaEoEua':
            re_string1 = re.compile('(Ea|Eo|Eu)')
            re_string2 = re.compile('a')
            return re_string2.sub('',re_string1.sub('E',p0)) == re_string2.sub('',re_string1.sub('E',p1))
        else:
            return False
    else:
        return False

'''Report statistics for relationship between corpus and dictionary'''
dictionary_stems, corpus_stems, common_stems = dictionary_sets(stem_dictionary, parsed_dictionary)

def make_full_analyses(dictionary=stem_dictionary, limit=100000):

    analyses = {}

    #for (count, raw_stem) in enumerate(common_stems):
    for (count, raw_stem) in enumerate(dictionary):
        if count > limit:
            break
        else:
            #print count, raw_stem
            analyses[raw_stem]={'parsings':[], 'solutions':set([])}
            for (prefix,body) in geminate_tildes(substitute_initial(parse_prefix(raw_stem))):
                for (root,pattern) in analyse(body):
                    if root_dictionary.has_key(root):
                        rootcount = root_dictionary[root]
                    else:
                        rootcount = 0
                    if pattern_dictionary.has_key(prefix+pattern):
                        patterncount = pattern_dictionary[prefix+pattern]
                    else:
                        patterncount = 0
                    if rootcount*patterncount==0:
                        comment = ''
                    else:
                        comment = '*'
                        analyses[raw_stem]['solutions'].add((root,prefix+pattern,10*np.log2(rootcount*patterncount)))
                    analyses[raw_stem]['parsings'].append((prefix,body,root,rootcount,prefix+pattern,patterncount,comment))
                    
        #print_entry(logfile, raw_stem, analyses[raw_stem]['parsings'])
    return analyses

def report_detailed_analyses(dictionary=stem_dictionary, limit=100000):

    analyses = {}

    #for (count, raw_stem) in enumerate(common_stems):
    for (count, raw_stem) in enumerate(dictionary):
        if count > limit:
            break
        else:
            #print count, raw_stem
            analyses[raw_stem]={'parsings':[], 'solutions':set([])}
            for (prefix,body,code) in geminate_tildes(substitute_initial(parse_prefix(raw_stem,''))):
                for (root,pattern,code) in analyse(body,code):
                    if root_dictionary.has_key(root):
                        rootcount = root_dictionary[root]
                    else:
                        rootcount = 0
                    if pattern_dictionary.has_key(prefix+pattern):
                        patterncount = pattern_dictionary[prefix+pattern]
                    else:
                        patterncount = 0
                    if rootcount*patterncount==0:
                        comment = ' '
                    else:
                        comment = '*'
                        analyses[raw_stem]['solutions'].add((root,prefix+pattern,10*np.log2(rootcount*patterncount),code))
                    analyses[raw_stem]['parsings'].append((prefix,body,root,rootcount,prefix+pattern,patterncount,comment,code))
                    
        print_entry(logfile, raw_stem, analyses[raw_stem]['parsings'])
        
        dictionarydic = make_parse_dictionary()
        if raw_stem in common_stems:
            (dic_root, dic_pattern) = dictionarydic[raw_stem]
            dic_root = str(dic_root)
            dic_pattern = str(dic_pattern)
            print >> logfile,'Dic      '+raw_stem+':     ('+dic_root+', '+dic_pattern+')\n'
        else:
            print >> logfile,'No dictionary entry for ',raw_stem,'\n'

    return analyses

def make_terse_analyses(dictionary=stem_dictionary, limit=100000):

    analyses = {}

    #for (count, raw_stem) in enumerate(common_stems):
    for (count, raw_stem) in enumerate(dictionary):
        if count > limit:
            break
        else:
            #print count, raw_stem
            analyses[raw_stem]={'parsings':[], 'solutions':set([])}
            for (prefix,body) in geminate_tildes(substitute_initial(parse_prefix(raw_stem))):
                for (root,pattern) in analyse(body):
                    if root_dictionary.has_key(root):
                        rootcount = root_dictionary[root]
                    else:
                        rootcount = 0
                    if pattern_dictionary.has_key(prefix+pattern):
                        patterncount = pattern_dictionary[prefix+pattern]
                    else:
                        patterncount = 0
                    if rootcount*patterncount==0:
                        comment = ''
                    else:
                        comment = '*'
                        analyses[raw_stem]['solutions'].add((root,prefix+pattern,10*np.log2(rootcount*patterncount)))
                    #analyses[raw_stem]['parsings'].append((prefix,body,root,rootcount,prefix+pattern,patterncount,comment))
                    
        #print_entry(logfile, raw_stem, analyses[raw_stem]['parsings'])
    return analyses

logfile = open(baampath+'BAAM-09-02-10-solutions.log', 'wt')

#analyses = make_full_analyses(dictionary=common_stems, limit=len(common_stems))
#testlist = ['muriyH','>us~','sum~','Tuwl','$uwraY']
testlist = common_stems
analyses = report_detailed_analyses(dictionary=testlist, limit = 250 )

#analyses = report_detailed_analyses(dictionary=['muriyH'])

##solutions_by_counts = {}
##
##for key in analyses.keys():
##    try:
##        number_of_solutions = len(analyses[key]['solutions'])
##    except:
##        number_of_solutions = 0
##    try:
##        solutions_by_counts[number_of_solutions].append(key)
##    except:
##        solutions_by_counts[number_of_solutions] = [key]
##
##for count in solutions_by_counts.keys():
##    print >> logfile, count, len(solutions_by_counts[count])
##
##dictionarydic = make_parse_dictionary()
##
##re_string1 = re.compile('(Ea|Eo|Eu)')
##re_string2 = re.compile('(a)')
##    
##total_listed = 0
##total_unlisted = 0
##total_agree = 0
##total_disagree = 0
##
##for count in solutions_by_counts.keys():
##
##    listed = 0
##    unlisted = 0
##    agree = 0
##    disagree = 0
##
##    if count == 0:
##        print >> logfile, 'No solutions found ...'
##        for unsolved in solutions_by_counts[0]:
##            if unsolved in common_stems:
##                (dic_root, dic_pattern) = dictionarydic[unsolved]
##                dic_root = str(dic_root)
##                dic_pattern = str(dic_pattern)
##                print >> logfile,unsolved ,': (', dic_root,', ',dic_pattern,')'
##                #print_entry(logfile, unsolved, analyses[unsolved]['parsings'])
##                listed += 1
##            else:
##                unlisted += 1
##    elif count == 1:
##        print >> logfile, 'One solution found ...'
##        for singleton in solutions_by_counts[1]:
##            (root,pattern,logfreq) = list(analyses[singleton]['solutions'])[0]
##            if singleton in common_stems:
##                (dic_root, dic_pattern) = dictionarydic[singleton]
##                dic_root = str(dic_root)
##                dic_pattern = str(dic_pattern)
##                #print (root,pattern,dic_root, dic_pattern,root==dic_root,re_string1.sub('E',pattern)==re_string1.sub('E',dic_pattern),root_pattern_match(root,pattern,dic_root, dic_pattern,strictness='EaEo'))
##                if root_pattern_match(root,pattern,dic_root, dic_pattern,strictness='EaEoEua'):
##                    agree += 1
##                else:
##                    print >> logfile, '[', singleton, '~ (', root,', ', pattern,') # (', dic_root,', ', dic_pattern,')]'
##                    #print_entry(logfile, singleton, analyses[singleton]['parsings'])
##                    disagree += 1
##    else:
##        print >> logfile, 'Several solutions found ... [',count,']'
##        for multiple_solutions in solutions_by_counts[count]:
##            if multiple_solutions in common_stems:
##                solution_list = list(analyses[multiple_solutions]['solutions'])
##                where = np.argmax([l for (r,p,l) in solution_list])
##                (root,pattern,logfreq) = solution_list[where]
##                (dic_root, dic_pattern) = dictionarydic[multiple_solutions]
##                dic_root = str(dic_root)
##                dic_pattern = str(dic_pattern)
##                if root_pattern_match(root,pattern,dic_root, dic_pattern,strictness='EaEoEu'):
##                    agree += 1
##                else:
##                    print >> logfile, '[', multiple_solutions, '~ (', root,', ', pattern,') # (', dic_root,', ', dic_pattern,')]'
##                    #print_solutions(logfile, multiple_solutions, analyses[multiple_solutions]['solutions'])
##                    disagree += 1
##    if count>0:
##        print (count,len(solutions_by_counts[count]), agree, disagree,(agree*100.)/(agree+disagree))
##    if count==0:
##        print (count,len(solutions_by_counts[count]), listed, unlisted)
##
##    total_listed += listed
##    total_unlisted += unlisted
##    total_agree  +=  agree
##    total_disagree += disagree
##
##
##print >> logfile, total_listed, 'Listed, ', total_unlisted, ' Unlisted'
##print >> logfile, total_agree, 'Agree, ', total_disagree, ' Disagree'
##print >> logfile, (total_agree*100.)/(total_agree+total_disagree), '% agreement'
## 
###anafile = open('analyses.pkl','wb')
###cPickle.dump(analyses,anafile)
##
###solfile = open('solutions_by_counts.pkl','wb')
###cPickle.dump(solutions_by_counts,solfile)
##
###anafile.close()
###solfile.close()
        
logfile.close()

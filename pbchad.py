##### Licensed under the GNU General Public License (Version 3.0)
##### See LICENSE for more details.  
##### Written by Alexander W. Neal Riasanovsky.  

from __future__ import print_function

from src import g_to_g6
#import networkx as nx
from numpy import matrix

import sys
import fileinput

if len(sys.argv) <= 2:
    K = 32
    SOLVER_NUM = 0

if len(sys.argv) == 3:
    K = int(sys.argv[2])
    SOLVER_NUM = 0

if len(sys.argv) >= 4:
    K = int(sys.argv[2])
    SOLVER_NUM = int(sys.argv[3])

if SOLVER_NUM == 0:
    from pysat.solvers import Glucose3
if SOLVER_NUM == 1:
    from pysat.solvers import Glucose4
if SOLVER_NUM == 2:
    from pysat.solvers import Lingeling
if SOLVER_NUM == 3:
    from pysat.solvers import MapleChrono
if SOLVER_NUM == 4:
    from pysat.solvers import MapleCM
if SOLVER_NUM == 5:
    from pysat.solvers import Maplesat
if SOLVER_NUM == 6:
    from pysat.solvers import Minicard
if SOLVER_NUM == 7:
    from pysat.solvers import Minisat22
if SOLVER_NUM == 8:
    from pysat.solvers import MinisatGH

from pysat.pb import *

DEBUGGING = False
VERBOSE = True

if DEBUGGING:
    import datetime
    import time


TO_HADS = "./hads/"
TO_CNFS = "./cnfs/"

def diagonal_matrix(x):
    m = len(x)
    d = []
    for i in range(m):
        currow = [0 for j in range(m)]
        currow[i] = x[i]
        d.append(currow)
    
    return matrix(d)



#from Steve
def swap_cols(M,i,j):
    m = M.shape[0]
    n = M.shape[1]      #M.dimensions()
    N = M.copy()
    x = [M[k,i]  for k in range(m)]        #M.column(i)
    for k in range(m):
        N[k,i]=N[k,j]
        N[k,j]=x[k]
    return N


#checks to see if c1 = A|B and c2 = A|~B|C somehow
#if not, returns -2; if b is empty formula, return -1; else, return index of B in c1
def simplify_helper(c1, c2):
    
    k = 0
    b = -1
    
    #search the elements of c1 one-by-one
    while k < len(c1):
        
        #if k is not a term of A...
        if not c1[k] in c2:
            
            #c2 can only miss one element of c1, namely B
            if b > -1 or not -c1[k] in c2:
                return -2
            
            #otherwise, define b
            b = k
        k += 1
    return b


#helper method for abc_simplifier
def abc_helper(c1, c2):
    k = 0
    b = -1
    while k < len(c1):
        
        #if k is not a term of A...
        if not c1[k] in c2:
            
            #c2 can only miss one element of c1, namely B
            if b > -1 or not -c1[k] in c2:
                return -2
            
            #otherwise, identify B
            b = k
        k += 1
    return b
        

#a more dynamic version of simplify_helper
#takes in the entire clause dictionary and mutates it if possible
def abc_simplifier(good_clauses):
    for i in good_clauses:
        c1 = good_clauses[i]
        for j in good_clauses:
            if i != j:
                c2 = good_clauses[j]
                helper_output = abc_helper(c1, c2)
                if helper_output == -1:
                    
                    #then c1 = A and c2 = A|C
                    #remove c2
                    
                    if DEBUGGING:
                        print('at', (i,j), c1, c2, '    are of the form A and A|C', file = sys.stderr)
                        print('so remove', c2, file = sys.stderr)
                    dud = good_clauses.pop(j)
                    abc_simplifier(good_clauses)
                    return True
                else:
                    if helper_output >= 0:
                        
                        #then c1 = A|B and c2 = A|~B|C
                        #replace c2 with A|C
                        
                        if DEBUGGING:
                            print('at', (i,j), c1, c2, '    are of the form A|B and A|~B|C', file = sys.stderr)
                            print('so remove', -c1[helper_output], 'from', c2, file = sys.stderr)
                        
                        c2.remove(-c1[helper_output])
                        good_clauses[j] = c2
                        abc_simplifier(good_clauses)
                        return True
    return False



#helper mutator method which relabels one variable at a time
#removes trivially satisfied clauses if spotted; 
#also removes dupes within clauses if spotted
def relabel(good_clauses, a, b):
    if DEBUGGING:
        print('    okay, relabeling', b, 'to', a, file = sys.stderr)
    to_pop = []
    for i in good_clauses:
        c1 = good_clauses[i]
        c1.sort(cmp = lambda x,y:abs(x)-abs(y))
        
        if a in c1:
            if -a in c1 or -b in c1:
                #clause is always satisfied
                to_pop.append(i)
                if DEBUGGING:
                    print('        c1=', c1, 'so toss it!', file = sys.stderr)
            else:
                if b in c1:
                    if DEBUGGING:
                        print('        c1=', c1, 'so remove', b, file = sys.stderr)
                    #removing b is an equivalent clause
                    c1.remove(b)
                    good_clauses[i] = c1
                #else, don't change anything
        else:
            if -a in c1:
                if b in c1:
                    if DEBUGGING:
                        print('        c1=', c1, 'so toss it!', file = sys.stderr)
                    #clause is always satisfied
                    to_pop.append(i)
                else:
                    if -b in c1:
                        if DEBUGGING:
                            print('        c1=', c1, 'so remove', -b, file = sys.stderr)
                        #removing -b is an equivalent clause
                        c1.remove(-b)
                        good_clauses[i] = c1
                    #else, don't change anything
            else:
                if b in c1:
                    if -b in c1:
                        if DEBUGGING:
                            print('        c1=', c1, 'so toss it!', file = sys.stderr)
                        #clause is always satisfied
                        to_pop.append(i)
                    else:
                        if DEBUGGING:
                            print('        c1=', c1, 'so swap', b, 'with', a, file = sys.stderr)
                        c1.remove(b)
                        c1.append(a)
                        c1.sort(cmp = lambda x,y: abs(x)-abs(y))
                        good_clauses[i] = c1
                else:
                    if -b in c1:
                        if DEBUGGING:
                            print('        c1=', c1, 'so swap', -b, 'with', -a, file = sys.stderr)
                        c1.remove(-b)
                        c1.append(-a)
                        c1.sort(cmp = lambda x,y:abs(x)-abs(y))
                        good_clauses[i] = c1
    for i in to_pop:
        dud = good_clauses.pop(i)
                        
            

#searches for pairs of clauses c1, c2 where:
    #c1 = x_a | ~x_b and c2 = ~x_a | x_b (a < b)
    #replaces x_b with x_a everywhere
    #if x_b is a dummy variable, remove c1 and c2
def simplify_by_substitutions(good_clauses, dummy_threshold):
    for i in good_clauses:
        c1 = good_clauses[i]
        if len(c1) == 2:
            for j in good_clauses:
                if i != j:
                    c2 = good_clauses[j]
                    if len(c2) == 2:
                        if c1[0]==-c2[0] and c1[1]==-c2[1]:
                            a = c1[0]
                            b = -c1[1]
                            if DEBUGGING:
                                print('at', (i,j),':', c1, c2, file = sys.stderr)
                            if abs(b) >= dummy_threshold:
                                dud = good_clauses.pop(i)
                                dud = good_clauses.pop(j)
                                if DEBUGGING:
                                    print('we dont need either clause now!', file = sys.stderr)
                            else:
                                #good_clauses[i] = c1
                                #good_clauses[j] = c2
                                if DEBUGGING:
                                    print('we need to keep the clauses intact', file = sys.stderr)
                            relabel(good_clauses, a, b)
                            simplify_by_substitutions(good_clauses, dummy_threshold)
                            return True
    return False


#helper method for finding gaps in a list
#requires that L is sorted
def has_gaps(L):
    for i in range(len(L)-1):
        if L[i+1] > L[i]+1:
            return i
    return -1

#shifts dummy variables if there are gaps
def dummy_shift(good_clauses, dummy_threshold):
    shifted = False
    #maxvar = 0
    dums = []
    for i in good_clauses:
        for k in good_clauses[i]:
            #if abs(k) > maxvar:
            #    maxvar = abs(k)
            if not abs(k) in dums:
                if abs(k) >= dummy_threshold:
                    dums.append(abs(k))
    dums.sort()
    i = has_gaps(dums)
    while i >= 0:
        shifted = True
        if DEBUGGING:
            print('    found a gap:', (dums[i], dums[i+1]), 'in', dums, file = sys.stderr)
        relabel(good_clauses, dums[i]+1, dums[i+1])
        dums[i+1] = dums[i] + 1
        i = has_gaps(dums)
    return shifted
        

def remove_singletons(good_clauses, dummy_threshold):
    for i in good_clauses:
        if len(good_clauses[i]) == 1:
            to_pop = []
            a = good_clauses[i][0]
            for j in good_clauses:
                if i != j:
                    if a in good_clauses[j]:
                        to_pop.append(j)
                        if DEBUGGING:
                            print(good_clauses[j], 'got popped by', good_clauses[i], file = sys.stderr)
                    if -a in good_clauses[j]:
                        if DEBUGGING:
                            print(good_clauses[j], 'got shortened by', good_clauses[i], file = sys.stderr)
                        new_clause = good_clauses[j]
                        new_clause.remove(-a)
                        good_clauses[j] = new_clause
            if abs(a) >= dummy_threshold:
                to_pop.append(i)
                if DEBUGGING:
                    print('at', i, good_clauses[i], 'is a dummy singleton', file = sys.stderr)
            for j in to_pop:
                dud = good_clauses.pop(j)
                print('    removed', a,'at', j, 'wow!', file = sys.stderr)
            remove_singletons(good_clauses, dummy_threshold)
            return True
    return False



#takes in a list of tuples corresponding to a CNF formula 
#returns a slightly shorter list for an equivalent formula
#if there exist clauses c1 c2 where: 
    #c1 = A|B and c2 = A|~B|C, 
    #then c1&c2 = (A|B)&(A|~B|C) <==> (A|B)&(A|C),
    #and if B or C is the empty formula, then c1&c2 <==> A
    #or if instead A is the empty formula, then c1&c2 <==> C

#if there exist two variables x_i, x_j (i<j) where: 
    #clauses x_i|~x_j and ~x_i|x_j exist, then replace x_j with x_i 
    #if j >= dummy_threshold, remove x_j entirely

#if there exists a singleton clause x_i (or ~x_i), then:
    #for all clauses c containing x_i, remove the entire clause
    #for all clauses c containing ~x_i, remove ~x_i from the clause

#once no such simplifications exist, shift all dummy variables 
#down to the lowest possible names for the dummy variables

#assumes no tuple in old_tups has both i and -i
def simplify_clause_list(old_clauses, dummy_threshold):
    
    #keep a dictionary of all good clauses, by their original index
    good_clauses = dict()
    if DEBUGGING:
        print('original clauses:', file = sys.stderr)
    #make sure the clauses are sorted
    for i in range(len(old_clauses)):
        curr_clause = old_clauses[i]
        curr_clause.sort(cmp = lambda x,y: abs(x)-abs(y))
        good_clauses[i] = curr_clause
        if DEBUGGING:
            print(i, curr_clause, file = sys.stderr)
    
    lctr = 0
    progressing = True
    while progressing:
        progressing = False
        
        if DEBUGGING:
            print('printing nonzeroed clauses...', file = sys.stderr)
            show_zeroed_clauses( good_clauses, [-i for i in range(1,dummy_threshold)] )
        
        
        removed_singles = remove_singletons(good_clauses, dummy_threshold)
        if DEBUGGING and removed_singles:
            print('removed singles!', file = sys.stderr)
        
            print('now only', len(good_clauses), 'clauses!', file = sys.stderr)
            print('printing nonzeroed clauses...', file = sys.stderr)
            show_zeroed_clauses( good_clauses, [-i for i in range(1,dummy_threshold)] )
        
        subbed = simplify_by_substitutions(good_clauses, dummy_threshold)
        if DEBUGGING and subbed:
            print('subbed!', file = sys.stderr)
        
            print('now only', len(good_clauses), 'clauses!', file = sys.stderr)
            print('printing nonzeroed clauses...', file = sys.stderr)
            show_zeroed_clauses( good_clauses, [-i for i in range(1,dummy_threshold)] )
        
        abc_simplified = abc_simplifier(good_clauses)
        if DEBUGGING and abc_simplified:
            print('abc simplified!', file = sys.stderr)

            print('now only', len(good_clauses), 'clauses!', file = sys.stderr)
            print('printing nonzeroed clauses...', file = sys.stderr)
            show_zeroed_clauses( good_clauses, [-i for i in range(1,dummy_threshold)] )
        
        progressing = subbed or abc_simplified or removed_singles
        lctr += 1 
        if DEBUGGING:
            print('finished', lctr, 'loops!', file = sys.stderr)
        
    shifted = dummy_shift(good_clauses, dummy_threshold)
    if DEBUGGING and shifted:
        print('shifted!', file = sys.stderr)
    
        print('now only', len(good_clauses), 'clauses!', file = sys.stderr)
        print('printing nonzeroed clauses...', file = sys.stderr)
        show_zeroed_clauses( good_clauses, [-i for i in range(1,dummy_threshold)] )
    
    maxvar = 0
    for i in good_clauses:
        for k in good_clauses[i]:
            if abs(k) > maxvar:
                maxvar = abs(k)
    return [good_clauses[i] for i in good_clauses], maxvar

#helper method for debugging; helps show if a clause is 
#satisfied by all-false truth assignments
def show_zeroed_clauses(good_clauses, assumed):
    print('assumed', assumed, file = sys.stderr)
    for i in good_clauses:
        reduced = [j for j in good_clauses[i] if not -j in assumed]
        j = 0
        zeroed = False
        
        if len(reduced) == 0:
            print('contradiction at', i, good_clauses[i], 'with assumption:', assumed, file = sys.stderr)
            return False
        
        while j < len(assumed) and not zeroed:
            if assumed[j] in reduced:
                zeroed = True
            j += 1
        
        if not zeroed:
            if len(reduced) == 1:
                print('   added', reduced[0], 'because of', good_clauses[i], file = sys.stderr)
                assumed.append(reduced[0])
                show_zeroed_clauses(good_clauses, assumed)
                return
    
    satisfied = True
    for i in good_clauses:
        reduced = [j for j in good_clauses[i] if not -j in assumed]
        j = 0
        zeroed = False
        while j < len(assumed) and not zeroed:
            if assumed[j] in reduced:
                zeroed = True
            j += 1
        
        if not zeroed and len(reduced) > 0:
            if DEBUGGING:
                print(i, reduced, file = sys.stderr)
            satisfied = False
    print('satisfied by zero?', satisfied, file = sys.stderr)
    print("", file = sys.stderr)
    return



#assumes K <= 62
def sol_to_g6(sol):
    curr_g6 = str(chr(K+63))
    curr_int = 0
    curr_index = 5
    
    
    
    is_bad = False
    nums = []
    correct_length = 1
    
    if not len(curr_g6) == correct_length:
        print('not the right length immediately!', curr_g6, file = sys.stderr)
    
    
    for k in sol:
        if k > 0:
            curr_int += (1<<curr_index)
        if curr_index == 0:
            curr_g6 += str(chr(curr_int+63))
            
            
            
            
            correct_length += 1
            if not len(curr_g6) == correct_length:
                print('not the right length!', curr_g6, file = sys.stderr)
    
            
            
            if curr_int+63 > 126:
                is_bad = True
            nums.append(curr_int+63)
            
            
            
            curr_index = 5
            curr_int = 0
        else:
            curr_index -= 1
    if curr_index < 5:
        curr_g6 += str(chr(curr_int+63))
        
        
        
        correct_length += 1
        if not len(curr_g6) == correct_length:
            print('not the right length!', curr_g6, file = sys.stderr)
        
        if curr_int+63 > 126:
            is_bad = True
        nums.append(curr_int+63)
        #print('lengthened to', curr_g6, 'and done', file = sys.stderr)
    
    
    if is_bad:
        print('bad!', nums, file = sys.stderr)
    
    
    if not len(curr_g6) == 84:
        print('not the right length!', curr_g6, file = sys.stderr)
    
    return curr_g6
    



#takes in a list of Hadamard matrices and returns a list of CNF tuples
#the CNF tuples encode the diagonalization problem as a SAT problem
def hads_to_graphs(all_columns = True, transpose = True):
    
    if DEBUGGING:
        start_time = datetime.datetime.now()
    
    translate = {'0':[-1,-1,-1,-1], '1':[-1,-1,-1, 1], '2':[-1,-1, 1,-1], '3':[-1,-1, 1, 1], '4':[-1, 1,-1,-1], '5':[-1, 1,-1, 1], '6':[-1, 1, 1,-1], '7':[-1, 1, 1, 1], '8':[ 1,-1,-1,-1], '9':[ 1,-1,-1, 1], 'A':[ 1,-1, 1,-1], 'B':[ 1,-1, 1, 1], 'C':[ 1, 1,-1,-1], 'D':[ 1, 1,-1, 1], 'E':[ 1, 1, 1,-1], 'F':[ 1, 1, 1, 1]}
    
    #tracks the number of graphs and hadamard matrices
    hctr = 0
    gctr = 0
    
    mults_to_clauses = dict()
    
    emap = dict()
    ectr = 0
    for a in range(K):
        for b in range(a+1,K):
            ectr += 1
            emap[ectr] = (a,b)
            emap[(a,b)] = ectr
    
    re_emap = dict()
    ectr = 0
    for b in range(K):
        for a in range(b):
            ectr += 1
            #re_emap[(a,b)] = ectr
            re_emap[ emap[(a,b)] ] = ectr
            print('relabeled index of', (a,b), 'from', emap[(a,b)], 'to', ectr, file = sys.stderr)
            
    #return
    
    if all_columns:
        max_col = K-1
        
    else:
        max_col = 0

    for curr_line in fileinput.input():
        myh = []
        #print('curr line is:')
        #print(curr_line)
        #print('just did it')
        #print('currline[:-1] is:')
        #print(curr_line[:-1])
        #print('just did it')
        #print()
        if DEBUGGING:
            print('the entries of curr line are:', file = sys.stderr)
            for s in curr_line:
                print(s)
        
        for s in curr_line[:-1]:
            if DEBUGGING:
                print(s, file = sys.stderr)
                print(translate[s], file = sys.stderr)
            myh.extend(translate[s])
        myh = matrix([myh[i*K:(i+1)*K] for i in range(K)])
            
        
        #for line in chunk.splitlines():
        #    myh = []
        #    for s in line:
        #        myh.extend(translate[s])
        #myh = matrix([myh[i*K:(i+1)*K] for i in range(K)]])
        
        if VERBOSE:
            print('trying matrix', hctr, '...', file = sys.stderr)
        
        for col in range(max_col+1):
            hh = swap_cols(myh, 0, col)
            dd = diagonal_matrix( [hh[j, 0] for j in range(K)])
            
            hh = dd*hh
            
            dd = diagonal_matrix( [hh[0,j] for j in range(K)] )
            hh = hh*dd

            if transpose:
                newhh = []
                for i in range(K):
                    newhh.append([hh[j,i] for j in range(K)])
                hh = matrix(newhh)
            
            nclauses = 0
            nvars = K*(K-1)/2
            
            if SOLVER_NUM == 0:
                g = Glucose3()
            if SOLVER_NUM == 1:
                g = Glucose4()
            if SOLVER_NUM == 2:
                g = Lingeling()
            if SOLVER_NUM == 3:
                g = MapleChrono()
            if SOLVER_NUM == 4:
                g = MapleCM()
            if SOLVER_NUM == 5:
                g = Maplesat()
            if SOLVER_NUM == 6:
                g = Minicard()
            if SOLVER_NUM == 7:
                g = Minisat22()
            if SOLVER_NUM == 8:
                g = MinisatGH()
            
            matchings = dict()
            
            ectr = 0
            
            for b in range(K-1):
                ectr += 1
                my_row = [0 for i in range(K-1)]
                my_row[b] = K
                my_row = tuple(my_row)
                matchings[my_row] = [(ectr,0,b)]
            
            for a in range(1,K):
                for b in range(a+1,K):
                    if DEBUGGING:
                        print('trying', (a,b), file = sys.stderr)
                    ectr += 1
                    hab =tuple([hh[a,j] * hh[b,j] for j in range(K)])
                    my_row = []
                    for j in range(1,K):
                        dot = 0
                        for k in range(K):
                            dot += hab[k]*hh[j,k]
                        my_row.append(dot)
                    
                    my_row = tuple(my_row)
                    
                    if my_row in matchings:
                        matchings[my_row].append((ectr,a,b))
                        old_ectr = matchings[my_row][0][0]
                        if DEBUGGING:
                            print('the old ectr,ectr=', old_ectr, ectr, file = sys.stderr)
                        g.add_clause((re_emap[old_ectr],-re_emap[ectr]))
                        g.add_clause((-re_emap[old_ectr],re_emap[ectr]))
                        nclauses += 2
                    
                    else:
                        matchings[my_row] = [(ectr,a,b)]
                        
                        coes_to_inds = dict()
                        mults = []
                        
                        
                        for i in range(K-1):
                            coe = my_row[i]
                            if coe in coes_to_inds:
                                coes_to_inds[coe].append(i)
                            else:
                                coes_to_inds[coe] = [i]
                        
                        mults = [ (coe,len(coes_to_inds[coe])) for coe in coes_to_inds if coe]
                        mults.sort()
                        mults = tuple(mults)
                        
                        if mults in mults_to_clauses:
                            old_ectr, num_new_vars, old_coes_to_inds, old_clauses = mults_to_clauses[mults]
                            
                            #must relabel the old indices to match the new ones
                            #for each coefficient c 
                            relab = dict()
                            
                            for coe in old_coes_to_inds:
                                coe_ctr = 0
                                for old_index in old_coes_to_inds[coe]:
                                    relab[1+old_index] = re_emap[1+coes_to_inds[coe][coe_ctr]]
                                    relab[-1-old_index] = -re_emap[1+coes_to_inds[coe][coe_ctr]]
                                    coe_ctr += 1
                            
                            
                            relab[old_ectr] = re_emap[ectr]
                            relab[-old_ectr] = -re_emap[ectr]
                            
                            for i in range(1,1+num_new_vars):
                                relab[old_ectr+i] = i+nvars
                                relab[-old_ectr-i] = -i-nvars
                            
                            for curr_clause in old_clauses:
                                g.add_clause([relab[k] for k in curr_clause])
                            
                            nclauses += len(old_clauses)
                            nvars += num_new_vars
                            
                        else:
                            if DEBUGGING:
                                print('dealing with', mults,'for the first time...', file = sys.stderr)
                            cnf = PBEnc.equals( lits = range(1,K) + [ectr], weights = list(my_row)+[-K], bound = 0, encoding = 4 )
                            curr_clauses, maxvar = simplify_clause_list(cnf.clauses, ectr + 1)
                            
                            num_new_vars = maxvar - ectr
                            
                            if DEBUGGING:
                                print("", file = sys.stderr)
                                print("", file = sys.stderr)
                                print('adding new clauses; there are', nvars, 'vars', file = sys.stderr)
                            
                            for curr_clause in curr_clauses:
                                if DEBUGGING:
                                    print('relabeling clause', c, 'with ectr=', ectr,'and', nvars, 'vars', file = sys.stderr)
                                new_clause = []
                                for k in curr_clause:
                                    if abs(k) <= ectr:
                                        if k > 0:
                                            new_clause.append(re_emap[k])
                                        else:
                                            new_clause.append(-re_emap[-k])
                                    else:
                                        if k > 0:
                                            if DEBUGGING:
                                                print(k, 'is the', k-ectr,'th dummy...', 'shift it up by', nvars, file = sys.stderr)
                                            new_clause.append(k-ectr+nvars)
                                        else:
                                            new_clause.append(k+ectr-nvars)
                                            if DEBUGGING:
                                                print(k, 'is the', k+ectr,'th dummy...', 'shift it down by', nvars, file = sys.stderr)
                                if DEBUGGING:
                                    print('    adding', new_clause, file = sys.stderr)
                                g.add_clause(new_clause)
                            
                            nvars += num_new_vars
                            nclauses += len(curr_clauses)
                            
                            mults_to_clauses[mults] = (ectr, num_new_vars, coes_to_inds, curr_clauses)
                            
                            if DEBUGGING:
                                print('now there are', (nvars, nclauses), 'vars,clauses', file=sys.stderr)
                            
            sol_ctr = 0
            
            if VERBOSE:
                print('there are', nvars, 'variables and', nclauses, 'clauses', file=sys.stderr)
                print('trying to find all solutions!',file=sys.stderr)
            
            if True:    #nvars >= 2000:
                while g.solve():
                    new_sol = g.get_model()
                    sol_ctr += 1
                    
                    if sol_ctr % 1000 == 0:
                        print('   ', sol_ctr, 'found so far...', file = sys.stderr)
                    
                    gctr += 1
                    #nx.write_graph6(nx.Graph([emap[k] for k in new_sol[:K*(K-1)/2] if k > 0]), sys.stdout, nodes = range(K))
                    #print( str([emap[k] for k in new_sol[:K*(K-1)/2] if k > 0]), file = sys.stdout )
                    g.add_clause( [-new_sol[1+j*(j+1)/2] for j in range(K-1)] )
                    #print([k for k in new_sol[:K*(K-1)/2] if k > 0], file = sys.stderr)
                    #print('which is really:', file = sys.stderr)
                    #print([re_emap[emap[k]] for k in new_sol[:K*(K-1)/2] if k > 0], file = sys.stderr)
                    #print('got written as\n', file = sys.stderr)
                    print(sol_to_g6(new_sol[:K*(K-1)/2]), file = sys.stdout)
                    #g.add_clause( [-new_sol[j] for j in range(K-1)] )
                
                if VERBOSE:
                    print(sol_ctr, 'solutions found!', file = sys.stderr)
                    if DEBUGGING:
                        end_time = datetime.datetime.now()
                        elapsed_time = end_time - start_time
                        print('total time elapsed:', elapsed_time.seconds,":",elapsed_time.microseconds, 'matrices solved:', hctr, file = sys.stderr)
            
            else:
                if VERBOSE:
                    print('too hard for now!', file = sys.stderr)
            
            hctr += 1
    return

hads_to_graphs(all_columns = False, transpose = True)

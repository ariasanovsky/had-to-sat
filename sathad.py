from itertools import combinations    
from src import *
from numpy import matrix


#debug variables; delete once code is good!
VERBOSE = True
VERBOSE_SIMP = False
VERBOSE_HADS = False
VERBOSE_MULTS = False
VERBOSE_LIST = False
VERBOSE_OR = False
VERBOSE_SOLS = True
TIDYING = False

TO_HADS = "./hads/"
TO_CNFS = "./cnfs/"
TO_CERTS = "./certs/"

TO_SLOANE = TO_HADS + "sloane/"
TO_KTR = TO_HADS + "ktr/"

def diagonal_matrix(x):
    m = len(x)
    d = []
    for i in range(m):
        currow = [0 for j in range(m)]
        currow[i] = x[i]
        d.append(currow)
    
    #print "diag of", x, "is:"
    #for r in d:
        #print r, 'of length', len(r)
    #print
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


#gets all auxiliary solutions corresponding to a multiset of coefficients
#for example if mults = ((-4, 3), (4, 6)) is the multiset 
#and the dimension K = 12, then -4*alpha+4*beta in {0, 12} 
#s.t. 0 <= alpha <= 3 and 0 <= beta <= 6 is the auxilairy equation
#we incldue the final coordinate indicating if it attained 0 or K
def all_sols(mults, K):
    all_ieqs = []
    final_ieq = []
    for i in range(len(mults)):#[0. .len(mults)-1]:
        coe,mul = mults[i]
        new_eq = [0 for j in range(len(mults))]#[0. .len(mults)]]
        new_eq[i] = 1
        all_ieqs.append( [0] + new_eq )
        new_eq[i] = -1
        all_ieqs.append( [mul] + new_eq )
        final_ieq.append( coe )
        
    new_eq = [0 for j in range(len(mults)+1)]#[0. .len(mults)]]
    new_eq[-1] = 1
    all_ieqs.append( [0] + new_eq )
    new_eq[-1] = -1
    all_ieqs.append( [1] + new_eq )
        
    all_ieqs.append([0] + [cc for cc in final_ieq] + [-K])
    all_ieqs.append([0] + [-cc for cc in final_ieq] + [K])
    
    if VERBOSE:
        print 'there are', len(all_ieqs), 'inequalities'
    pp = Polyhedron( ieqs = all_ieqs )
    sols = pp.integral_points()
    
    if VERBOSE:
        print 'YOU SHOULD NOT BE HERE'
        print 'there are', len(sols), 'solutions'
        print 'they are:' 
        for sol in sols:
            print '    ', sol
    return sols













def test_equivalence( tups1, tups2, num_vars, num_tries = 1000 ):
    ctr = 0
    subgen = Subsets(range(1,num_vars+1))  #[1. .num_vars])
    while ctr < num_tries:
        ctr += 1
        
        sub = subgen.random_element()
        compl = Set( [j for j in range(1,num_vars+1) if not j in sub]) #[1. .num_vars] if not j in sub] )
        
        t1_satisfied = True
        i = 0
        while t1_satisfied and i < len(tups1):
            si_pos = Set([j for j in tups1[i] if j>0])
            si_neg = Set([abs(j) for j in tups1[i] if j<0])
            t1_satisfied = len(sub.intersection(si_pos))>0 or len(compl.intersection(si_neg))>0
            i+=1
        
        t2_satisfied = True
        i = 0
        while t2_satisfied and i < len(tups2):
            si_pos = Set([j for j in tups2[i] if j>0])
            si_neg = Set([abs(j) for j in tups2[i] if j<0])
            t2_satisfied = len(sub.intersection(si_pos))>0 or len(compl.intersection(si_neg))>0
            i+=1
        
        if t1_satisfied != t2_satisfied:
            return False
    return True








#takes in a list of tuples corresponding to a CNF formula 
#returns a shorter list for an equivalent formula
#only simplifies if it finds clauses c1 c2 where: 
    #c1 = a|b and c2 = a|~b|c
    #generally, c1&c2 = (a|b)&(a|~b|c) <==> (a|b)&(a|c)
    #if b or c is null, then c1&c2 <==> a

#it may assumed that the CNF formula has no tautologicaly clauses
#i.e., no tuple in old_tups has both i and -i
def simplify_tup_list(old_tups):
    
    irredundants = dict()
    
    for i in range(len(old_tups)):      #[0. .len(old_tups)-1]:
        irredundants[i] = set(old_tups[i])
    
    progressing = True
    while progressing:
        if VERBOSE_SIMP:
            print 'num good rows:', len(old_tups)
        progressing = False
        for i in range(len(old_tups)):#[0. .len(old_tups)-1]:
            if i in irredundants:
                for j in range(i+1, len(old_tups)):#[i+1. .len(old_tups)-1]:
                    if j in irredundants and i in irredundants:
                        sj = irredundants[j]
                        si = irredundants[i]
                        inters = si.intersection(sj)
                        diffi = si.difference(sj)
                        diffj = sj.difference(si)
                        
                        
                        #first check if one is a subset of the other
                        #if so, delete the bigger one
                        if len(diffi) == 0:
                            if VERBOSE_SIMP:
                                print si, 'is a subset of', sj
                            dead_ind = irredundants.pop(j)
                            progressing = True
                        else:
                            if len(diffj) == 0:
                                if VERBOSE_SIMP:
                                    print sj, 'is a subset of', si
                                dead_ind = irredundants.pop(i)
                                progressing = True
                            else:
                                
                                
                                #checking to see if rows i and j 
                                #are of the form a|b, and a|~b|c
                                #where b is a single variable in diffi or diffj
                                #and a is the "intersection" of the clauses
                                #note that (a|b)&(a|~b|c) <==> (a|b)&(a|c)
                                #and if c is vacuous, then instead we have 
                                #(a|b)&(a|~b) <==> a
                                if len(diffi) == 1:
                                    b = list(diffi)[0]
                                    if -b in sj:
                                        if len(diffj) == 1:
                                            if VERBOSE_SIMP:
                                                print si, sj, 'are of the form (a|b)&(a|~b) <==> a'
                                                print 'here b=',b
                                            dead_index = irredundants.pop(j)
                                            si.remove(b)
                                            irredundants[i] = si
                                            progressing = True
                                        else:
                                            if VERBOSE_SIMP:
                                                print si, sj, 'are of the form (a|b)&(a|~b|c) <==> (a|b)&(a|c)'
                                                print 'here b=',b
                                            sj.remove(-b)
                                            irredundants[j] = sj
                                            progressing = True
                                else:
                                    if len(diffj) == 1:
                                        b = list(diffj)[0]
                                        if -b in si:
                                            if VERBOSE_SIMP:
                                                print sj, si, 'are of the form (a|b)&(a|~b|c) <==> (a|b)&(a|c)'
                                            si.remove(-b)
                                            irredundants[i] = si
                                            progressing = True
    new_tups = []
    for i in irredundants:
        new_tups.append( tuple(irredundants[i]) )
    return new_tups













#naive way of computing the Cartesian product of a list of sets
def cart_prod(S):
    if len(S) == 0:
        return []
    
    curr_prods = [[j] for j in S[0]]
    
    #print 'first set of trivial tuples:', curr_prods
    #print 
    i = 1
    while i < len(S):
        prev_prods = curr_prods
        
        #print 'previous tuples:', prev_prods
        #print 'merging with trivial tuples:', [[j] for j in S[i]]
        
        curr_prods = []
        for prev_prod in prev_prods:
            for j in S[i]:
                curr_prods.append( prev_prod+[j] )
        i+=1
    return [prev_prod for prev_prod in curr_prods]


#naive method for find all solutions to the auxiliary equations
#of the form -4*alpha + 4*beta + . . . = 0 OR K
def all_sols2(mults, K):
    sols = []
    coe_vec = tuple( [mult[0] for mult in mults] )
    poss_sols = cart_prod([range(mults[j][1]+1) for j in range(len(mults))])
    #cart_prod([[0. .mults[j][1]] for j in [0. .len(mults)-#1]])
    
    #if VERBOSE_SOLS:
    #    print 'there are', len(poss_sols), 'possible solutions'
    
    for poss_sol in poss_sols:
        dot = 0        #dot = coe_vec * vector(poss_sol)
        for j in range(len(coe_vec)):
            dot += coe_vec[j] * poss_sol[j]
        #print 'dot of', coe_vec,'and', vector(poss_sol),'is', dot
        if dot == 0:
            sols.append(tuple(poss_sol + [0]))
            #print 'added 0'
            #print
        if dot == K:
            sols.append(tuple(poss_sol + [1]))
            #print 'added K'
            #print
    #print 'but only', len(sols), 'solutions'
    if VERBOSE_SOLS:
        print '    the', len(sols), 'sols:'
        print sols
    
    return sols



















#efficiently OR's together a list of clauses
def or_recur( clauses ):
    cnum = len(clauses)
    if cnum:
        if cnum == 1:
            #if VERBOSE_OR:
                #print '    one clause, length', len(str(clauses[0]))
            return [clauses[0]]
        
        if cnum == 2:
            #if VERBOSE_OR:
                #print 'two clauses, lengths', len(str(clauses[0])), len(str(clauses[1]))
            c0 = clauses[0] | clauses[1]
            
            
            #no idea which cnf converter is most efficient
            
            if VERBOSE_OR:
                print 'before convert_cnf_recur', len(str(c0))
            
            #c0.convert_cnf()
            c0.convert_cnf_recur()
            
            t0 = CNF_to_tup_list(c0)
                        
            t1 = simplify_tup_list(t0)
            c1 = tup_list_to_CNF(t1)
            
            if VERBOSE_OR:
                print len(t0), 'clauses to', len(t1)
                print len(str(c0)), 'characters to', len(str(c1))
            
            #thesat = c0.satformat()
            
            #if VERBOSE_OR:
                #print 'this should be in CNF already... should be trivial, no?'
            
            #if VERBOSE_OR:
                #print 'shouldnt change the num clauses...', len(t2)
                #print 'length is similar?', len(str(c2))
            
            
            return [c1]
        
        c0s = clauses[:cnum/2]
        c1s = clauses[cnum/2:]
        
        return or_recur(or_recur(c0s) + or_recur(c1s))
    return []


#takes the aux solutions and multiset and writes it as a CNF formula
#where each variable is a dummy which counts the number of 1s next to 
#a particular coefficient
def sols_to_CNF(sols, mults):
    
    #first, we will make len(sols) AND clauses, 
    #one for each solution corresponding to mults
    sol_clauses = []
    
    num_coes = len(mults)
    
    dummap = dict()
    dummy_ctr = 1
    
    for i in range(num_coes):#[0. .num_coes-1]:
        dummap[i] = [j+dummy_ctr for j in range(mults[i][1]+1)]#[0. .mults[i][1]] ]
        dummy_ctr += mults[i][1]+1
    
    #if VERBOSE_SOLS:
    #    print "the dummies are", dummap
    
    #sols = list(sols)
    #sols.sort()
    
    for sol in sols:
        #print sol, dummap, mults
        gg = propcalc.formula('y'+str(dummap[0][sol[0]]))
        for i in range(1,num_coes):#[1. .num_coes-1]:
            idums = dummap[i]
            gg = gg & propcalc.formula('y'+str(idums[sol[i]]))
            
            #gg.convert_cnf_recur()
        
        #if VERBOSE_SOLS:
            #print 'solution', sol, 'corresponds to', gg
        sol_clauses.append(gg)
    
    sol_clause = or_recur(sol_clauses)[0]
    
    #sol_clause.convert_cnf()
    
    #print 'simplifying sol_clause of length', len(str(sol_clause))
    #sol_clause = tup_list_to_CNF(simplify_tup_list(CNF_to_tup_list(sol_clause)))
    #print 'go sol_clause of length', len(str(sol_clause))
    
    #we need clauses to force exactly one dummy variable 
    #for each coefficient to be true
    
    
    
    #do dum_clause in hads_to_tup_list instead! 
    
    
    
    #dum_clause = or_recur([propcalc.formula('y'+str(dum)) for dum in dummap[0]])[0]
    
    
    #if VERBOSE:
        #print 'first dummy:', dum_clause
    
    #for j in [0. .mults[0][1]]:
        #for k in [j+1. .mults[0][1]]:
            #curr = propcalc.formula("~y"+str(dummap[0][j])+"|~y"+str(dummap[0][k]))
            #if VERBOSE:
                #print 'pair dummy:', curr
            #dum_clause = dum_clause & curr
            #dum_clause.convert_cnf_recur()



    #for i in [1. .num_coes-1]:
        #idums = dummap[i]
        #curr = or_recur([propcalc.formula('y'+str(dum)) for dum in idums])[0]
        #dum_clause = dum_clause & curr
        
        #if VERBOSE:
            #print 'next dummy:', curr
        #dum_clause.convert_cnf_recur()
        
        #for j in [0. .mults[i][1]]:
            #for k in [j+1. .mults[i][1]]:
                #curr = propcalc.formula("~y"+str(idums[j])+"|~y"+str(idums[k]))
                #if VERBOSE:
                    #print 'pair dummy:', curr
                #dum_clause = dum_clause & curr
                #dum_clause.convert_cnf_recur()
    
    #test_sols(sols, dum_clause & sol_clause, dummap)
    
    #return [sol_clause, dum_clause, dummap]
    return [sol_clause, dummap]
    

#takes a multiset of coefficients and writes it as a list of tuples
#where each tuple is a row of the associated SAT file, in terms of 
#the dummy variables which count ones in front of a coefficient
def mults_to_tup_list(mults, K):
    tups = []
    [sol_clause, dummap] = sols_to_CNF(all_sols2(mults, K), mults)
    #print "WOW GOT HERE"
    
    for clause in str(sol_clause).split("&"):
        #print '    clause:', clause, 'is a', type(clause)
        tup = []
        #print '    clause[1:-1]:', clause[1:-1], 'is a', type(clause[1:-1])
        #print '    str(clause[1:-1]):', str(clause[1:-1]), 'is a', type(str(clause[1:-1]))
        
        if VERBOSE_MULTS:
            print 'the clause is', clause
        
        for va in str(clause[1:-1]).split("|"):
            if va[0] == 'y':
                tup.append(int(va[1:]))
            else:
                tup.append( -int(va[2:]) )
        
        if VERBOSE_MULTS:
            print 'the tup is', tup
            print 
        
        tups.append(tuple(tup))
    
    #for clause in str(dum_clause).split("&"):
        #print '    clause:', clause, 'is a', type(clause)
        #tup = []
        #print '    clause[1:-1]:', clause[1:-1], 'is a', type(clause[1:-1])
        #print '    str(clause[1:-1]):', str(clause[1:-1]), 'is a', type(str(clause[1:-1]))
        #for va in str(clause[1:-1]).split("|"):
            #if va[0] == 'y':
                #tup.append(int(va[1:]))
            #else:
                #tup.append( -int(va[2:]) )
        #tups.append(tuple(tup))
    
    return [tups, dummap]


#takes in a list of tuples and writes a CNF file
def tup_list_to_SAT( tup_list, nvars ):
    the_SAT = str('p cnf '+str(nvars) + ' ' + str(len(tup_list)) + '\n')
    for tup in tup_list:
        curr = str()
        for num in tup:
            curr = curr + str(num) + " "
        the_SAT = the_SAT + curr[:-1] + " 0 \n"
    return the_SAT




def tup_list_to_CNF(tup_list):
    
    form = str()
    
    for tup in tup_list:
        curr = '('

        for i in range(len(tup)):
            if tup[i]>0:
                curr = curr + 'y' + str(tup[i]) + '|'
            else:
                curr = curr + '~y' + str(abs(tup[i])) + '|'
        if VERBOSE_LIST:
            print 'curr is', curr
        form = form + curr[:-1] + ')&'
    if VERBOSE_LIST:
        print 'form looks like:'
        print form[:-1]
        print
    
    form = propcalc.formula(form[:-1])
    if VERBOSE_LIST:
        print 'now form looks like:'
        print form
        print
    
    return form




def CNF_to_tup_list(form):
    tups = []
    for clause in str(form).split("&"):
        #print '    clause:', clause, 'is a', type(clause)
        tup = []
        #print '    clause[1:-1]:', clause[1:-1], 'is a', type(clause[1:-1])
        #print '    str(clause[1:-1]):', str(clause[1:-1]), 'is a', type(str(clause[1:-1]))
        
        #if VERBOSE_LIST:
            #print 'the clause is', clause
        
        for va in str(clause[1:-1]).split("|"):
            if va[0] == 'y':
                tup.append(int(va[1:]))
            else:
                tup.append( -int(va[2:]) )
        
        if VERBOSE_LIST:
            print 'the tup is', tup
            print 
        
        tups.append(tuple(tup))
    return tups






#takes in a list of Hadamard matrices and returns a list of CNF tuples
#the CNF tuples encode the diagonalization problem as a SAT problem
def hads_to_tups(hads, all_columns = True, get_cheats = False, get_matrices = False, get_tups = False, had_type = 'sloane'):
    
    
    #entries L_H are lists corresponding to each hadamard matrix H
    #L_H contains L_{H, col}, one list for each choice of the normalized column col 
    #L_{H, c} contains tuples which will be a line of the CNF file
    
    
    if get_cheats:
        cheatsheets = []
    
    if get_matrices:
        mats = []
    
    if get_tups:
        master_tup_lists = []
    
    
    
    
    
    
    if had_type == 'sloane':
        sgn_to_num = lambda x: 1 if x=='+' else -1
    
    if had_type == 'ktr':
        translate = {'0':[-1,-1,-1,-1], '1':[-1,-1,-1, 1], '2':[-1,-1, 1,-1], '3':[-1,-1, 1, 1], '4':[-1, 1,-1,-1], '5':[-1, 1,-1, 1], '6':[-1, 1, 1,-1], '7':[-1, 1, 1, 1], '8':[ 1,-1,-1,-1], '9':[ 1,-1,-1, 1], 'A':[ 1,-1, 1,-1], 'B':[ 1,-1, 1, 1], 'C':[ 1, 1,-1,-1], 'D':[ 1, 1,-1, 1], 'E':[ 1, 1, 1,-1], 'F':[ 1, 1, 1, 1]}
        myhad = hads[0]
        infile = open(TO_KTR + myhad, "r")
        infile.readline()
        infile.readline()
        infile.readline()
        infile.readline()
        nchars = 288 #32 rows, 8 hex chars, one \n
        chunk = infile.read(nchars)
        infile.readline()
        infile.readline()
        infile.readline()
    
    #master list which takes a row constraint multiset and 
    #stores the tup list associated with it, in terms of dummy variables
    tup_map = dict()
    h_ctr = 0
    
    c = 0
    
    while c < len(hads):           #for myhad in hads:
        if had_type == 'sloane':
            myhad = hads[c]
            infile = open(TO_SLOANE + myhad, "r") 
            myh = []
            for line in infile:
                if len(line)>0 and line[0] in ['+', '-']:
                    #print 'the line:', line
                    myh.append([sgn_to_num(ent) for ent in line if ent in ['+', '-']])
            infile.close()
            myh = matrix(myh)
        
        
        if had_type == 'ktr':
            myh = []
            if not chunk:
                myhad = hads[c]
                infile = open(TO_KTR + myhad, "r")
                infile.readline()
                infile.readline()
                infile.readline()
                infile.readline()
                chunk = infile.read(nchars)
                infile.readline()
                infile.readline()
                infile.readline()
            
            for line in chunk.splitlines():
                currline = []
                for s in line:
                    currline.extend(translate[s])
                    #print 'extending by', translate[s], '...'
                myh.append(currline)
            #print 'ended with length', len(myh)
            
            myh = matrix(myh)
        
        if had_type == 'sage':
            myh = hadamard_matrix_www(myhad)
        
        
        K = myh.shape[0]
        max_col = K-1
        
        if not all_columns:
            max_col = 0
        
        tup_lists = []
        
        if get_cheats:
            curr_cheats = []
        
        if get_matrices:
            curr_mats = []
        
        
        if VERBOSE:
            print 'trying matrix', myhad, '. . .'
        
        for col in range(max_col+1):#[0. .max_col]:
            if VERBOSE and all_columns:
                print 'trying column', col, '. . .'
            
            hh = swap_cols(myh, 0, col)
            
            
            #dd = diagonal_matrix( hh.column(0))
            
            
            dd = diagonal_matrix( [hh[j, 0] for j in range(K)])
            
            hh = dd*hh
            #print
            #print "after first diagonalization:"
            #print hh
            
            dd = diagonal_matrix( [hh[0,j] for j in range(K)] )
            hh = hh*dd
            #print "after second diagonalization:"
            #print hh
            #print
            
            matchings = dict()
            
            curr_clauses = []
            
            num_clauses = 0
            
            if get_cheats:
                cheatsheet = dict()
                dum_cheat = dict()
            
            if get_matrices:
                C = []
            
            #in the end, the SAT file has 1, 2, . . ., K-1 representing 
            #l_{12}, . . ., l_{1,K} and the rest are dummy variables
            #there is a set of dummies for each row of C
            #on each row, there is a subset of variables for each coefficient c
            #write m for the multiplicity of c
            #these m+1 dummies indicate the number of 1's in front of c
            
            if get_cheats:
                for i in range(1,K):#[1 K-1]:
                    cheatsheet[i] = "x" + str(i) + ": edge 1," + str(i+1)
            
            dummy_ctr = K-1
            
            #print hh
            #for a in range(K): #[0. .K-1]:
                #print [hh[a]*hh[b] for b in #range(K)]     #[0. .K-1]]
                #print vector( [hh[a,j] * hh[0,j] for j in range(K)])       #[0. .K-1]] ) == hh[a]
            
            
            
            outname = "./cnfs/h"+str(K)+"n" + str(h_ctr) + "c" + str(col) + ".cnf"

            #print outname
            
            outfile = open(outname,"w")
            outfile.close()
            
            
            
            for a in range(K):      #[0. .K-1]:
                for b in range(a+1,K):   #[a+1. .K-1]:
                    hab =tuple([hh[a,j] * hh[b,j] for j in range(K)])
                           #vector( [hh[a,j] * hh[b,j] for j in range(K)])    #[0. .K-1]] )
                    #print 'for a,b=', (a,b),'hab=', hab
                                #my_row = tuple([hab*hh[j] for j in range(1,K)])     #[1. .K-1]])
                    my_row = []
                    for j in range(1,K):
                        dot = 0
                        for k in range(K):
                            dot += hab[k]*hh[j,k]
                        my_row.append(dot)
                    my_row = tuple(my_row)
                    
                    #if VERBOSE:
                        #print my_row
                    if my_row in matchings:
                        matchings[my_row].append((a,b))
                    else:
                        matchings[my_row] = [(a,b)]
                        
                        if get_matrices:
                            C.append(my_row)
                        
                        #print 'new row at:', (a,b),":", my_row
                        
                        if a > 0:
                            #print '    processing row:', (a,b), '. . .'
                            coes = []
                            inds = dict()
                            for j in range(K-1):        #[0. .K-2]:
                                coe = my_row[j]
                                if coe:
                                    if coe in inds:
                                        inds[coe].append(j)
                                    else:
                                        inds[coe] = [j]
                                        coes.append(coe)
                            coes.sort()
                            mults = tuple( [(coe, len(inds[coe])) for coe in coes] )
                            
                            if not mults in tup_map:
                                if True:#VERBOSE_HADS:
                                    print '    making a new tup list out of', mults
                                    
                                tup_map[mults] = mults_to_tup_list(mults, K)
                                if True:#VERBOSE_HADS:
                                    print '    added a tup map of', len(tup_map[mults][0])
                                    #print
                                
                                #if VERBOSE_HADS:
                                    #print '        now the known mults are', list(tup_map)
                                    #print '        made the tup list!'
                                    #print 
                            
                                #if VERBOSE_HADS:
                                    #print '        the known mults are', list(tup_map)
                                
                                
                                #if VERBOSE_HADS:
                                    #print '        the known mults are', list(tup_map)

                            [row_clauses, dummap] = tup_map[mults]
                            
                            dum_relabel = dict()
                            new_dumbs = 0
                            for i in dummap:
                                if VERBOSE_HADS:
                                    print 'relabeling y', dummap[i], 'to y', [j+dummy_ctr for j in dummap[i]]
                                
                                coe = mults[i][0]
                                for j in range(len(dummap[i])):     #[0. .len(dummap[i])-1]:
                                    yind = dummap[i][j]
                                    
                                    dum_relabel[yind] = yind+dummy_ctr
                                    dum_relabel[-yind] = -yind-dummy_ctr
                                    
                                    
                                    if get_cheats:
                                        cheat = 'y_(' + str(a) + "," + str(b) + "); (" + str(coe) + "," + str(j) + "):"
                                        cheat = cheat + " num of " + str(coe) + " xs which are "
                                        cheat = cheat + str(["x"+str(k+1) for k in inds[coe]]) + " =" + str(j)
                                    
                                        cheatsheet[yind+dummy_ctr] = cheat
                                    
                                        dum_cheat[yind+dummy_ctr] = (j, [k+1 for k in inds[coe]])
                                
                                
                                new_dumbs += len(dummap[i])
                            
                            dummy_ctr += new_dumbs
                            
                            
                            
                            with open(outname,"a") as outfile:
                                for tup in row_clauses:
                                    new_line = tuple([dum_relabel[j] for j in tup])
                                    curr_clauses.append(new_line)
                                    num_clauses+=1
                                    strs=" ".join(str(x) for x in new_line)
                                    outfile.write(strs+" 0\n")

                            
                            
                            #pretty sure the issue is in the cardinality constraints. . .
                            #we tested the solution clause already. . .
                            
                            if VERBOSE_HADS:
                                print 
                                print 'all inds are', inds
                            
                            with open(outname, "a") as outfile:
                                for i in range(len(mults)):     #[0. .len(mults)-1]:
                                    idums = dummap[i]
                                    coe, mult = mults[i]
                                    
                                    if VERBOSE_HADS:
                                        print 'coe,mult', coe, mult
                                        print 'indices', inds[coe]
                                        print 'dummies', idums
                                        print
                                    
                                    
                                    #at least one y_j is true
                                    #####doing the dummy_clause stuff here now
                                    
                                    if not len(idums) == mult+1:
                                        print 'index disagreement!'
                                        print 'indices', idums
                                        print 'coe,mult', coe, mult
                                        print 
                                    new_line = tuple([dum_relabel[idums[j]] for j in range(mult+1)])
                                    curr_clauses.append(new_line)      #[0. .mult]]))
                                    num_clauses+=1
                                    
                                    strs=" ".join(str(x) for x in new_line)
                                    outfile.write(strs+" 0\n")

                                
                                
                                
                                    for j in range(mult+1):     #[0. .mult]:
                                        
                                        
                                        #no two y_j's can be simultaneously true
                                        ####doing the dummy_clause stuff here now
                                        
                                        for k in range(j+1, mult+1):#[j+1. .mult]:
                                            curr_clauses.append( (-dum_relabel[idums[j]], -dum_relabel[idums[k]]) )
                                            num_clauses+=1
                                            outfile.write(str(-dum_relabel[idums[j]]) + " " + str(-dum_relabel[idums[k]]) + " 0\n")
                                        
                                        #y_j implies that every j+1 set of x's has a False
                                        for sub in combinations( inds[coe], j+1):
                                            tup = tuple([-w-1 for w in sub] + [dum_relabel[-idums[j]]])
                                            if VERBOSE_HADS:
                                                print '    ', tup
                                            curr_clauses.append(tup)
                                            num_clauses+=1
                                            strs=" ".join(str(x) for x in tup)
                                            outfile.write( strs+" 0\n" )
                                    
                                    
                                        #y_j implies that every mult-j+1 set of x's has a True
                                        for sub in combinations( inds[coe], mult-j+1):
                                            tup = tuple([w+1 for w in sub] + [dum_relabel[-idums[j]]])
                                            if VERBOSE_HADS:
                                                print '    ', tup
                                            curr_clauses.append(tup)
                                            num_clauses+=1
                                            strs=" ".join(str(x) for x in tup)
                                            outfile.write( strs+" 0\n" )

                                
            if TIDYING:# or len(curr_clauses)<11000:
                if VERBOSE:
                    print 'final tidying up of the', len(curr_clauses),'clauses. . .'
                simpler = simplify_tup_list(curr_clauses)
                tup_lists.append(simpler)
            
                if VERBOSE:
                    print '    added tup list of', len(simpler), 'constraint and', dummy_ctr, 'vars from', len(matchings), 'matchings'
                    print 
            else:
                tup_lists.append(curr_clauses)
                if VERBOSE:
                    print '    added tup list of', len(curr_clauses), 'constraints and', dummy_ctr, 'vars from', len(matchings), 'matchings'
                    print 
            
            if get_cheats:
                curr_cheats.append((cheatsheet, dum_cheat))
            
            if get_matrices:
                curr_mats.append(C)
            
        if get_tups:
            master_tup_lists.append(tup_lists)
        
                
        with open(outname, "r+") as outfile:
            existing=outfile.read()
            outfile.seek(0) #point to first line
            outfile.write("p cnf " + str(dummy_ctr) + " " + str(len(curr_clauses)) + "\n"+existing)
        
        
        if get_cheats:
            cheatsheets.append(curr_cheats)
        
        if get_matrices:
            mats.append(curr_mats)
        
        h_ctr += 1
        if not had_type == 'ktr':
            c += 1
        else:
            chunk = infile.read(nchars)
            if chunk:
                infile.readline()
                infile.readline()
                infile.readline()
            else:
                c += 1
                infile.close()
        
    if get_cheats:
        if get_matrices:
            if get_tups:
                return (master_tup_lists, mats, cheatsheets)
            else:
                return (mats, cheatsheets)

        else:
            if get_tups:
                return (master_tup_lists, mats)
            else:
                return mats
    
    else:
        if get_matrices:
            if get_tups:
                return (master_tup_lists, mats)
            else:
                return mats
        else:
            if get_tups:
                return master_tup_lists
            else:
                return
            
            
            
#outs = hads_to_tups(["./hads/had.16.0.txt", "./hads/had.16.1.txt", "./hads/had.16.2.txt", "./hads/had.16.3.txt", "./hads/had.16.4.txt"], all_columns = False)

#outs = hads_to_tups(["./hads/had.20.pal.txt", "./hads/had.20.will.txt", "./hads/had.20.toncheviv.txt"], all_columns = False)

#outs = hads_to_tups(["./hads/had.24." + str(j+1) + ".txt" for j in range(60)], all_columns = False)

#outs = hads_to_tups(["./hads/had.28." + str(j+1) + ".txt" for j in range(487)], all_columns = False)

#outs = hads_to_tups(["had.32." + hname + ".txt" for hname in  ["pal","syl","t1","t2","t3","t4"] ], all_columns = False)

outs = hads_to_tups(["ktr_test.txt" ], all_columns = False, had_type = 'ktr')


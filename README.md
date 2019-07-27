# had-to-sat
Licensed under the GNU General Public License (Version 3.0).  
See LICENSE for more details.  


Used to find all graphs diagonalizable by given Hadamard matrices.  
First converts the problem for each given Hadamard matrix to a SAT formula where edges (1,2),...,(K-1,K) correspond to variables x_1,...,x_{K(K-1)/2}.  
Each pair (a,b) (coresponding to x_j) of vertices yields a pseudo-boolean constraint of the form a_1x_1 + ... + a_{K-1}x_{K-1} = x_j.  
The CNF formula for pair (a,b) is then computed, simplified, and memoized according to the multiset of coefficients {a_1,...,a_{K-1}}.  
Then uses a SAT solver to find all solutions and writes them as graph6 strings.  Third option from the terminal selects one of 9 possible SAT solvers.  
Finally, use nauty to remove isomorphs.  

To set up:  
    1. install networkx and pysat with pip install networkx and pip install python-sat  
    2. Populate a source directory /src/ with https://searchcode.com/codesearch/raw/89954425/ saved as g_to_g6.py  
    3. Populate a directory /hads/ktr/ with Hadamard matrices from http://math.ipm.ac.ir/~tayfeh-r/Hadamard32.htm   
    4. Install nauty in a subdirectory.  

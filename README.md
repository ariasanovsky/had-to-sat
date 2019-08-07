# had-to-sat
Licensed under the GNU General Public License (Version 3.0).  
See LICENSE for more details.  


Used to find all graphs diagonalizable by given Hadamard matrices.  
First converts the problem for each given Hadamard matrix to a SAT formula where edges (0,1), (0,2), (1,2), ..., (K-2,K-1) correspond to variables x_1,...,x_{K(K-1)/2}.  Note: this is backwards lexicographical order which is more compatible with graph6.  
Each pair (a,b) (coresponding to x_j) of vertices yields a pseudo-boolean constraint of the form a_1x_1 + ... + a_{K-1}x_{K-1} = x_j.  
The CNF formula for pair (a,b) is then computed, simplified, and memoized according to the multiset of coefficients {a_1,...,a_{K-1}}.  
Then uses a SAT solver to find all solutions and writes them as graph6 strings.  Optional parameters can be used to select a solver (by default, Glucose3) and to specify a number of vertices other than 32.  
Finally, use nauty to remove isomorphs.  

To set up:  
    1. install pysat with pip install python-sat and set up nauty with make and configure.  
    2. Download Hadamard matrices from http://math.ipm.ac.ir/~tayfeh-r/Hadamard32.htm to a subdirectory and unrar them.  Recommendation: omit the Sylvester matrix, H32typethree.txt.   
    3. Process the Hadamard matrices so that they each appear as a line of 256 hexadecimal characters.  
    4. Pass this file through pbchad.py with python pbchad.py filename
    5. Refine output with nauty command shortg.  

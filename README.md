# had-to-sat
Licensed under the GNU General Public License (Version 3.0).  
See LICENSE for more details.  


Used to find all graphs diagonalizable by given Hadamard matrices.  
First converts the problem for each given Hadamard matrix to a SAT formula and writes to a cnf file.  
Then uses (any) SAT solver to find all solutions and write the given solutions as adjacency matrices of graphs.  
Then (not yet implemented), uses nauty to list all diagonalized graphs.  

To set up:  
    1. Populate a source directory /src/ with the contents of /logic/ as well as /misc/flatten.py from the sage /src/ directory.  
    2. Changing line 147 of boolformula.py so that it reads: from . import flatten  
    3. Populate a directory /hads/ktr/ with Hadamard matrices from http://math.ipm.ac.ir/~tayfeh-r/Hadamard32.htm   
    4. Populate a directory /hads/sloane/ with Hadamard matrices from http://neilsloane.com/hadamard/  
    5. Install your favorite SAT solver in the project directory.  
    6. (not implemented) Install nauty in a subdirectory.  

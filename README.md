# had-to-sat
Licensed under the GNU General Public License (Version 3.0).  
See LICENSE for more details.  


Used to find all graphs diagonalizable by given Hadamard matrices.  
First converts the problem for each given Hadamard matrix to a SAT formula and writes to a cnf file.  
Then uses any SAT solver (currently Glucse 4.0) to find all solutions and write the given solutions as graph6 strings.  
Then uses nauty to remove isomorphs.  

To set up:  
    1. install networkx and pysat with pip install networkx and pip install python-sat
    2. Populate a source directory /src/ with https://searchcode.com/codesearch/raw/89954425/ saved as g_to_g6.py  
    3. Populate a directory /hads/ktr/ with Hadamard matrices from http://math.ipm.ac.ir/~tayfeh-r/Hadamard32.htm   
    4. Install nauty in a subdirectory.  

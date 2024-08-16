# SDP_bounds_for_quantum_codes

Denote by $(( n, K, \delta ))$ a qubit quantum code with distance $\delta$, which encodes a K-dimensional Hilbert space to a $n$ qubit Hilbert space. 
  
This library provides a solution with a positive objective function of the SDP in Eq. (148) for  $((7,1,4))$ and of the SDP in Eq. (150) for $((8, 9, 3))$ and $(( 10 , 5, 4 ))$.

The solutions in pickle format: *inf\_cert\_7*, *inf\_cert\_8* and *inf\_cert\_10*. 

The programs *certificate\_7.ipynb*, *certificate\_8.ipynb* and *certificate\_10.ipynb* read the pickle files and check the objective function and the violation constraints for each three codes respectivly.

This proof the non-existence of qubit quantum codes $((7,1,4))$, $((8, 9, 3))$ and $(( 10 , 5, 4 ))$.

We additionally provide *lovasz.ipynb* which checks that the quantum lovasz number written as in Eq. (157), disproof the existence of the qubit quantum code $(( 4 , 1, 3 ))$.

All needed functions are defined in *fun_dual.ipynb*.
  

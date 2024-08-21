# SDP_bounds_for_quantum_codes

This work complements the results from https://arxiv.org/abs/2408.10323. All references are from the mentioned manuscript.

Denote by $(( n, K, \delta ))$ a qubit quantum code with distance $\delta$, which encodes a K-dimensional Hilbert space to a $n$ qubit Hilbert space. 
  
This library provides a solution with a positive objective function of the SDP in Eq. (148) for  $((7,1,4))$ and of the SDP in Eq. (150) for $((8, 9, 3))$ and $(( 10 , 5, 4 ))$.

The solutions are in pickle format: *inf\_cert\_7*, *inf\_cert\_8* and *inf\_cert\_10*. 

The programs *certificate\_7.ipynb*, *certificate\_8.ipynb* and *certificate\_10.ipynb* read the pickle files and check the following parameters for each three codes: the objective value and the violation constraints as well as the minimum eigenvalue.

This proves the non-existence of qubit quantum codes $((7,1,4))$, $((8, 9, 3))$ and $(( 10 , 5, 4 ))$.

The *SDP_dual.ipynb* includes the used two SDP programs. Both are written using PICOS (https://picos-api.gitlab.io/picos/) and they can be tunned in order to offer alternative certificates for $((7,1,4))$, $((8, 9, 3))$ and $(( 10 , 5, 4 ))$.

We additionally provide *lovasz.ipynb* which runs the SDP of the quantum Lovász number for self-dual codes written in Eq. (70). For the case of the qubit quantum code $(( 4 , 1, 3 ))$, this program shows that the maximum value is $\vartheta(G')=7$ and thus, $1+\vartheta(G') < 2^7$ leading to the non-existence of such code [see Corollary 9].

All needed functions are defined in *fun_dual.ipynb* and *qubit_upper_bounds.py* contains the information of the qubit quantum code table in https://www.codetables.de.
  

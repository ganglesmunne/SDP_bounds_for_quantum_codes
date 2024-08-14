from qiskit.quantum_info import pauli_basis,Pauli
import numpy as np
from itertools import product,permutations
import mpmath as mp
from numpy.linalg import eig
from scipy.linalg import block_diag
from collections import defaultdict
import time

mp.prec = 1000
mp.dps = 1000

#### Operations to Paulis ####

def wt(P,n):
    return n-list(P).count('I')
def supp_int(P1,P2,n):

    P1=list(P1)
    P2=list(P2)

    supp1 = set([k for k in range(n) if P1[k] != 'I'])
    supp2 = set([k for k in range(n) if P2[k] != 'I'])

    return len(supp1.intersection(supp2))
def eq_nontr(P1,P2,n):
    P1 = list(P1)
    P2 = list(P2)

    eq = [k for k in range(n) if P1[k] == P2[k] and P1[k] != 'I']

    return len(eq)


#### Functions ####
def multinomial(n,i,j,k,l):

    # value=mp.factorial(n)/(mp.factorial(i)*mp.factorial(j)*mp.factorial(k)*mp.factorial(l)*mp.factorial(n-i-j-k-l))
    # return int(value)

    value=1/(mp.factorial(i)*mp.factorial(j)*mp.factorial(k)*mp.factorial(l)*mp.factorial(n-i-j-k-l))
    return float(value)
def comb(n, k):

    value = mp.binomial(n, k)
    return int(value)
def gamma(i,j,t,p,n):
    return 3**(i+j-t)*2**(t-p)*multinomial(n,p,t-p,i-t,j-t)
def beta(i,j,k,m,t):

	b=0

	for u in range(m+1):

		b+= comb(u,t)*comb(m-2*k,u-k)*comb(m-k-u,i-u)*comb(m-k-u,j-u)*(-1)**(t-u)

	return b
def alpha(i,j,t,p,a,k,n):

    sum=0

    for g in range(0,p+1):
        sum += (-1)**(a-g)*comb(a,g)*comb(t-a,p-g)*2**(t-a-p+g)

    return beta(i-a,j-a,k-a,n-a,t-a)*3**((1/2)*(i+j)-t)*sum
def kron_delta(a,b):
    if a == b:
        return 1
    else:
        return 0
def K(m,k,n,ga,de):
    sum=0
    for alp in range(n+1):
        sum += (-1)**alp*comb(n-k,m-alp)*comb(k,alp)*ga**((n-k)-(m-alp))*de**(m-alp)
    return sum

### dic ###
def f4_fun(n):
    funct={}
    size = 0
    for i in range(n + 1):
        for j in range(i + 1):
            for t in range((i + j - n), min(i, j) + 1):
                for p in range(t + 1):
                    funct.update({(i, j, t, p): size})
                    size += 1
    return funct, size
def f2M_fun(n):
    fun={}
    count = 0
    for a in range(n+1):
        for k in range(a,int((n+a)/2)+1):
            fun.update({(a, k): count})
            count += 1
    return fun, count

#### sym ####

def xabtp(G, n):

    paulis = pauli_basis(n, True)

    E= {ind:p.to_label() for ind, p in enumerate(paulis)}

    l4n = range(4**n)

    x=[]
    f4={}
    count=0

    #c_list=[]

    for a in range(n+1):
        for b in range(a+1):
            for t in range((a+b-n),min(a,b)+1):
                for p in range(t+1):
                    sum=0

                    for i,j in product(l4n,l4n):

                        E1 = E[i]
                        E2 = E[j]

                        if wt(E1,n) == a and wt(E2,n) == b and supp_int(E1,E2,n) == t and eq_nontr(E1,E2,n) == p:

                            sum += G[i,j]
                            #c_list.append((i,j))

                    x.append(sum/gamma(a,b,t,p,n))
                    f4.update({(a,b,t,p): count})

                    count += 1
    return x, f4, E
def gammaprime(x, E, f4, n):

    gammap=np.zeros((4**n,4**n))

    l4n = range(4**n)

    for a in range(n+1):
        for b in range(a+1):
            for t in range((a+b-n),min(a,b)+1):
                for p in range(t+1):

                    abtp = f4[(a, b, t, p)]

                    M=np.zeros((4**n,4**n))

                    for i in l4n:
                        for j in range(i+1):

                            E1 = E[i]
                            E2 = E[j]

                            if wt(E1,n) == a and wt(E2,n) == b and supp_int(E1,E2,n) == t and eq_nontr(E1,E2,n) == p:

                                M[i,j] = 1
                    gammap += x[abtp]*(M+M.T-np.diag(np.diag(M)))

    return gammap
def alphamatrices(x, f4, n):

    f2M = {}
    Mat=[]
    count = 0

    for a in range(n+1):
        for k in range(a,int((n+a)/2)+1):

            matrix = [[0] * (n+a-2*k+1) for _ in range(n+a-2*k+1)]

            for i in range(k, n+a-k+1):
                for j in range(k, i + 1):

                    sum = 0

                    for t in range((i+j-n),min(i,j) + 1):
                        for p in range(t+1):

                            norm=(comb(n+a-2*k,i-k)*comb(n+a-2*k,j-k))**(1/2)
                            sum += alpha(i,j,t,p,a,k,n)*x[f4[(i, j, t, p)]]/norm
                    matrix[i-k][j-k] = sum
                    matrix[j-k][i-k] = sum

            Mat.append(matrix)
            f2M.update({(a,k):count})
            count += 1
    return Mat, f2M
def xabtp_zeros(G, n):

    paulis = pauli_basis(n, True)

    E= {ind:p.to_label() for ind, p in enumerate(paulis)}

    l4n = range(4**n)
    f4={}
    x=[]
    count=0

    for a in range(n+1):
        for b in range(a+1):
            for t in range((a+b-n),min(a, b) + 1):
                for p in range(t + 1):

                    sum=0

                    if ((t-p) % 2) == 1:
                        x.append(0)
                    else:
                        for i,j in product(l4n,l4n):

                            E1 = E[i]
                            E2 = E[j]

                            if wt(E1,n) == a and wt(E2,n) == b and supp_int(E1,E2,n) == t and eq_nontr(E1,E2,n) == p:

                                sum += G[i,j]

                        x.append(sum/gamma(a,b,t,p,n))

                    f4.update({(a,b,t,p): count})
                    count +=1

    return x, f4, E


### Enumerators ###

def B(A,D,n):

    B=[]

    for i in range(n+1):
        sum=0
        for k in range(n+1):
            sum += K(i,k,n,1,D**2-1)*A[k]
        B.append(sum/D**n)

    return B
def S(A,D,n):

    S=[]

    for i in range(n+1):
        sum=0
        for k in range(n+1):
            sum += (-1)**k*K(i,k,n,D-1,D+1)*A[k]
        S.append(sum/D**n)

    return S

def QMDS_weights(n,d,D):
    k = n-2*(d-1)
    a = (n+k)/2

    A = np.zeros(n+1)
    for j in range(n+1):
        s = 0
        for i in range(j+1):
            s = s + comb(j,i)*(-1)**(j-i) * D**(i-min(2*a-i,i))
        A[j] = D**(2*k)*comb(n,j)* s
    return A,k

## Check test ##
def random_hermitian(n):
    G = np.random.rand(4 ** n, 4 ** n)
    return G + G.conj().T
def random_positive(n):
    G = np.random.rand(4 ** n, 4 ** n)
    return G@G.conj().T

def permutation_group(n,d=4):
    perm_n = list(permutations(range(n)))
    perm_d = list(permutations(range(d)))
    perm_0 = list(permutations(range(1,d)))
    perm_0 = [tuple([0]+list(p)) for p in perm_0]
    return perm_n, perm_d, perm_0

def fun_gamma_tilde(gamma,n,dic,dic_inv,E,d=4):

    perm_n, perm_d, perm_0=permutation_group(n,d)

    size=range(d**n)

    gamma_p=np.zeros((d**n,d**n))

    select_local_permutation=list(product(range(int(mp.factorial(d-1))), repeat=n))

    for a,b in product(size,size):

        Ea=dic[a]
        Eb=dic[b]

        for full_perm in perm_n:

            Ea_perm=tuple([Ea[i] for i in full_perm])
            Eb_perm=tuple([Eb[i] for i in full_perm])

            for select_local in select_local_permutation:

                Ea_perm2 = tuple([perm_0[select_local[i]][Ea_perm[i]] for i in range(n)])
                Eb_perm2 = tuple([perm_0[select_local[i]][Eb_perm[i]] for i in range(n)])

                gamma_p[dic_inv[Ea]][dic_inv[Eb]] += gamma[dic_inv[Ea_perm2]][dic_inv[Eb_perm2]]

    return gamma_p

def fun_gamma_bar(gamma,n,dic,dic_inv,E,d=4):

    perm_n, perm_d, perm_0=permutation_group(n,d)

    size=range(d**n)

    gamma_t=np.zeros((d**n,d**n))

    select_local_permutation=list(product(range(int(mp.factorial(d))), repeat=n))

    for a,b in product(size,size):

        Ea=dic[a]
        Eb=dic[b]

        for full_perm in perm_n:

            Ea_perm=tuple([Ea[i] for i in full_perm])
            Eb_perm=tuple([Eb[i] for i in full_perm])

            # Ea_perm,Eb_perm= [Ea[i],Eb[i] for i in full_perm]

            # Ea_perm,Eb_perm

            for select_local in select_local_permutation:

                Ea_perm2 = tuple([perm_d[select_local[i]][Ea_perm[i]] for i in range(n)])
                Eb_perm2 = tuple([perm_d[select_local[i]][Eb_perm[i]] for i in range(n)])

                gamma_t[dic_inv[Ea]][dic_inv[Eb]] += gamma[dic_inv[Ea_perm2]][dic_inv[Eb_perm2]]

                # a1=dic_inv[Ea]
                # b1=dic_inv[Eb]
                # a2=dic_inv[Ea_perm]
                # b2=dic_inv[Eb_perm]
                # if (wt(E[a1], n) != wt(E[a2], n)) and (wt(E[b1], n) != wt(E[b2], n)) and (supp_int(E[a1], E[b1], n) != supp_int(E[a2], E[b2], n))  and (eq_nontr(E[a1], E[b1], n)!= eq_nontr(E[a2], E[b2], n)):
                #
                #     print([wt(E[a1],n),wt(E[a2],n)],[wt(E[b1],n),wt(E[b2],n)],[supp_int(E[a1], E[b1], n),supp_int(E[a2], E[b2],n)],
                #               [eq_nontr(E[a1], E[b1], n),eq_nontr(E[a2], E[b2],n)])

    return gamma_t

def Pauli_to_string(e):

    e=list(e)

    st=[]

    for i in e:
        if i=='I':
            st.append(0)
        elif  i=='X':
            st.append(1)
        elif i == 'Y':
            st.append(2)
        elif i == 'Z':
            st.append(3)
    return tuple(st)

def fun_lambda_k(x,f4,n):

    lambda_k=[]

    for k in range(n+1):

        sum1=0

        for a in range(n+1):
            for b in range(a+1):
                for t in range((a+b-n),min(a,b)+1):
                    for p in range(t+1):

                        abtp = f4[(a, b, t, p)]

                        if k == a + b - t - p:

                            sum1 += (2-kron_delta(a,b))*x[abtp]*gamma(a,b,t,p,n)


        lambda_k.append(sum1)

    return lambda_k

def gammatilde(x, E, f4, n):

    lambda_k,gamma_k=fun_lambda_k(x, f4,n)

    gammat=np.zeros((4**n,4**n))

    l4n = range(4**n)

    for k in range(n+1):

        M=np.zeros((4**n,4**n))

        for i in l4n:
            for j in range(i+1):

                E1 = E[i]
                E2 = E[j]

                if wt(E1,n) + wt(E2,n)-supp_int(E1,E2,n)-eq_nontr(E1,E2,n) == k:
                    M[i,j] = 1

        #gammat += (lambda_k[k]/gamma_k[k])*(M+M.T-np.diag(np.diag(M)))
        gammat += (lambda_k[k] / (4**n*3**k*comb(n,k))) * (M + M.T - np.diag(np.diag(M)))
    return gammat

def alphamatrices_complementary(x, f4, n):

    f2M = {}
    Mat=[]
    count = 0

    for a in range(n+1):
        for k in range(a,int((n+a)/2)+1):

            matrix = [[0] * (n+a-2*k+1) for _ in range(n+a-2*k+1)]

            for i in range(k, n+a-k+1):
                for j in range(k, i + 1):

                    sum = 0

                    for t in range((i+j-n),min(i,j) + 1):
                        for p in range(t+1):

                            norm=(comb(n+a-2*k,i-k)*comb(n+a-2*k,j-k))**(1/2)

                            sum += alpha(i,j,t,p,a,k,n)*(x[f4[(i+j-t-p,0,0,0)]]-x[f4[(i, j, t, p)]])/norm

                    matrix[i-k][j-k] = sum
                    matrix[j-k][i-k] = sum

            Mat.append(matrix)
            f2M.update({(a,k):count})
            count += 1
    return Mat, f2M

def alphamatrices_full(x, f4, n, K):

    f2M = {}
    Mat=[]
    count = 0

    for a in range(n+1):
        for k in range(a,int((n+a)/2)+1):

            matrix = [[0] * (n+a-2*k+1) for _ in range(n+a-2*k+1)]

            for i in range(k, n+a-k+1):
                for j in range(k, i + 1):

                    sum = 0

                    for t in range((i+j-n),min(i,j) + 1):
                        for p in range(t+1):

                            norm=(comb(n+a-2*k,i-k)*comb(n+a-2*k,j-k))**(1/2)

                            sum += alpha(i,j,t,p,a,k,n)*(x[f4[(i+j-t-p,0,0,0)]])/norm

                    matrix[i-k][j-k] = sum
                    matrix[j-k][i-k] = sum

            Mat.append(matrix)
            f2M.update({(a,k):count})
            count += 1
    return Mat, f2M

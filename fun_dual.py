import mpmath as mp
from itertools import product
import picos as pic
from itertools import product, permutations
import numpy as np

# binomials and coefficients:
def comb(n, k):
    """
    returns the binomial coefficients
    """

    value = mp.binomial(n, k)
    return int(value)
def multinomial(n,i,j,k,l):
    """
    returns the multinomial coefficients defined in Eq.(77)
    """
    value=mp.factorial(n)/(mp.factorial(i)*mp.factorial(j)*mp.factorial(k)*mp.factorial(l)*mp.factorial(n-i-j-k-l))

    return int(value)
def gamma(i,j,t,p,n):
    """
    returns the gamma coefficients defined in Eq.(76)
    """

    return 3**(i+j-t)*2**(t-p)*multinomial(n,p,t-p,i-t,j-t)
def beta(i,j,k,m,t):
    """
    returns the beta coefficients defined in Eq.(83)
    """

    return sum([comb(u,t)*comb(m-2*k,u-k)*comb(m-k-u,i-u)*comb(m-k-u,j-u)*(-1)**(t-u) for u in range(m+1)])

def alpha(i,j,t,p,a,k,n):
    """
    returns the alpha coefficients defined in Eq.(83)
    """
    sum=0

    for g in range(p+1):
        sum += (-1)**(a-g)*comb(a,g)*comb(t-a,p-g)*2**(t-a-p+g)

    return beta(i-a,j-a,k-a,n-a,t-a)*3**((1/2)*(i+j)-t)*sum

# mappings of (i,j,t,p) and (a,k) to 1D array.

def f4_fun(n):
    """
    dictionary that maps a 4-array (i,j,t,p) to a 1-D array.
    It defines the set I(q,n) in Eq.(75)
    """

    funct={}
    size = 0
    for i in range(n + 1):
        for j in range(n + 1):
            for t in range((i + j - n), min(i, j) + 1):
                for p in range(t + 1):
                    funct.update({(i, j, t, p): size})
                    size += 1
    return funct, size
def f2M_fun(n):
    """
    dictionary that maps a 2-array (a,k) to a 1-D array.
    Numerates the matrices from the block-diagonalization
    of the Terwilliger algebra [see Eq.(84)]
    """

    fun={}
    count = 0
    for a in range(n+1):
        for k in range(a,int((n+a)/2)+1):
            fun.update({(a, k): count})
            count += 1
    return fun, count

### dual yijtp ###
def yijtp_fun(Y,n):

    """
    returns the dual variable yijtp defined after Eq.(151)
    """

    f2M,_=f2M_fun(n)

    yijtp=[]

    for i,j in product(range(n + 1),range(n + 1)):
        for t in range((i + j - n), min(i,j) + 1):
            for p in range(t + 1):
                y=0
                for k in range(min(i,j)+1):
                    for a in range((max(i,j)+k-n),k+1):
                        if a>=0:
                            y +=alpha(i,j,t,p,a,k,n)*Y[f2M[a,k]][i-k,j-k]
                yijtp.append(y/gamma(i, j, t, p, n))
    return yijtp

def wt(P,n):

    """
    returns the weight of a pauli string P, e.g. 'XIZI', of size n.
    """
    return n-list(P).count('I')

def K_fun(j,i,n):
    """
    returns the quarternary Krawtchouk Polynomial as defined in Eq.(16) .
    """
    sum=0
    for alp in range(n+1):
        sum += (-1)**alp*comb(n-i,j-alp)*comb(i,alp)*3**(j-alp)
    return sum

def D_fun(C,n,d):
    """
    returns the Di function as defined in Eq.(150)) .
    """
    return [2 ** (-n) * sum(K(j, i, n) * C[j] for j in range(d)) for i in range(n + 1)]

def SDPdual(n, K, d, upperbound=None, solver='cvxopt', verbosity=0, sol_primal=True, sol_dual=True, sdptol=1e-09):

    f4, size = f4_fun(n)
    f2M, _ = f2M_fun(n)

    sdp = pic.Problem(verbosity=verbosity)
    sdp.options["*_tol"]=sdptol

    Y = [pic.SymmetricVariable("Y%s%s" % (a, k), int(n + a - 2 * k + 1)) for a in range(n + 1) for k in
         range(a, int((n + a) / 2) + 1)]

    C = pic.RealVariable("C", n+1)
    Q = pic.RealVariable("P", n + 1)
    D = [pic.sum([2**(-n)*K_fun(j,i,n)*C[j] for j in range(d)]) for i in range(n+1)]

    y = yijtp_fun(Y, n)

    sdp.add_list_of_constraints(Y[i] >> 0 for i in range(len(Y)))

    obj_fun = (K*D[0]-C[0])+(2**n/K-1)*Q[0] - y[f4[0, 0, 0, 0]]

    sdp.set_objective('max', obj_fun)

    if upperbound != None:
        sdp.add_constraint(obj_fun <= upperbound)

    #### Constraints ####

    sdp.add_list_of_constraints(
        -2*y[f4[i, 0, 0, 0]]-y[f4[i, i, i, i]] + K*D[i]-C[i]+(2**n/K-2)*Q[i]-Q[0]==0
        for i in range(1,d))

    sdp.add_list_of_constraints(
        -2*y[f4[i, 0, 0, 0]]-y[f4[i, i, i, i]] + (K*D[i]+(2**n/K-2)*Q[i]-Q[0])==0
        for i in range(d,n+1))

    sdp.add_list_of_constraints(
        pic.sum([y[f4[a, b, c, q]] + Q[a+b-c-q] for a, b, c, q in f4.keys()
                 if t - p == c - q and (i, j, i + j - t - p) in list(permutations((a, b, a + b - c - q)))])==0
    for i, j, t, p in f4.keys() if (t - p) % 2 == 0 and not j == 0 and not i==j==t==p)

    sol = sdp.solve(solver=solver, primals=sol_primal, duals=sol_dual)

    Y = [np.array(Y[i].value) if np.array(Y[i].value).size != 1 else np.array([[Y[i].value]]) for i in
                range(len(Y))]

    C, Q = C.value, Q.value
    obj=sdp.objective.value
    sol= sol.problemStatus

    return sol, Y, C, Q, f4, obj
import mpmath as mp
from itertools import product

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
    returns the dual variable yijtp defined after Eq.(164) or Eq.(174)
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
                yijtp.append(y)
    return yijtp

def wt(P,n):

    """
    returns the weight of a pauli string P, e.g. 'XIZI', of size n.
    """
    return n-list(P).count('I')

def K(j,i,n):
    """
    returns the quarternary Krawtchouk Polynomial as defined in Eq.(16) .
    """
    sum=0
    for alp in range(n+1):
        sum += (-1)**alp*comb(n-i,j-alp)*comb(i,alp)*3**(j-alp)
    return sum

def D_fun(C,n,d):
    """
    returns the Di function as defined after Eq.(174)) .
    """
    return [2 ** (-n) * sum(K(j, i, n) * C[j] for j in range(d)) for i in range(n + 1)]

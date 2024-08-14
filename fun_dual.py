import mpmath as mp
from itertools import product

# binomials and coefficients:

def comb(n, k):

    value = mp.binomial(n, k)
    return int(value)
def multinomial(n,i,j,k,l):

    value=mp.factorial(n)/(mp.factorial(i)*mp.factorial(j)*mp.factorial(k)*mp.factorial(l)*mp.factorial(n-i-j-k-l))

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

    for g in range(p+1):
        sum += (-1)**(a-g)*comb(a,g)*comb(t-a,p-g)*2**(t-a-p+g)

    return beta(i-a,j-a,k-a,n-a,t-a)*3**((1/2)*(i+j)-t)*sum

# mappings of i,j,t,p to 1D array and a,k to 1D array
def f4_fun(n):
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
    fun={}
    count = 0
    for a in range(n+1):
        for k in range(a,int((n+a)/2)+1):
            fun.update({(a, k): count})
            count += 1
    return fun, count

### dual yijtp ###
def yijtp_fun(Y,n):

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

# weight
def wt(P,n):
    return n-list(P).count('I')

def K(m,k,n,ga,de):
    sum=0
    for alp in range(n+1):
        sum += (-1)**alp*comb(n-k,m-alp)*comb(k,alp)*ga**((n-k)-(m-alp))*de**(m-alp)
    return sum

def D_fun(C,n,d):
    return [2 ** (-n) * sum(K(j, i, n, 1, 3) * C[j] for j in range(d)) for i in range(n + 1)]

# import pickle
#
# n,d=7,4
# decimals=9
#
# with open("inf_cert", 'rb') as inf_cert:
#     Y=pickle.load(inf_cert)
#
# y=yijtp_fun(Y,n)
#
# f4,_=f4_fun(n)
#
# ### obj value ###
#
# c = 2*y[f4[(n, 0, 0, 0)]]/gamma(n, 0, 0, 0, n) + y[f4[n, n, n, n]]/gamma(n, n, n, n, n)
#
# obj_val = -y[f4[0, 0, 0, 0]]/gamma(0, 0, 0, 0, n) - (2**n-1)*c
#
# print("obj_val=%s" % obj_val)
#
# #### Lovasz constraints ####
#
#
# const = True
#
# for i, j, t, p in f4.keys():
#     if not ((t - p) % 2) == 1 and (i + j - t - p) >= d and not j == 0 and not i==0:
#         zero = np.round(y[f4[i, j, t, p]], decimals)
#         if zero != 0:
#             const = False
#             print((i, j, t, p),zero)
#
# print("Lovasz zeros = %s" % const)
#
# const = True
#
# for i in range(d, n + 1):
#     zero = 2*y[f4[i, 0, 0, 0]] / gamma(i, 0, 0, 0, n) + y[f4[i, i, i, i]] / gamma(i, i, i, i, n) - c
#     zero = np.round(zero, decimals)
#     if zero != 0:
#         const = False
#         print(zero)
#
# print("Constraint zeros = %s" % const)
#
# ### minimal eigenvalue ###
#
# ref = np.infty
# for i in range(len(Y)):
#     min_eig = min(np.linalg.eigvals(Y[i].astype('float')))
#     if ref >= min_eig:
#         ref = min_eig
#
# print("Min eigenvalue = %s" % ref)

s = SymmetricFunctions(ZZ).s()

#Complement of a partition in the k x (n-k) rectangle.
def complement(k,n,eta):
    eta = eta + [0]*(k-len(eta))
    return [n-k-p for p in eta[::-1]]

#Partitions with n-1 boxes in k < (n-k).
def partitions(k,n):
    return [p+[0]*(k-len(p)) for p in Partitions(n-1) if (len(p) <= k and p[0] <= n-k)]

#Ribbon Schur function of composition b.
def ribbonSchur(b):
    k = len(b)
    bsum = [sum(b[:s]) for s in range(1,k+1)]
    lam = [bsum[k-i-1] - k + i + 1 for i in range(k)]
    mu = [bsum[k-i-1] - k + i for i in range(1,k)] + [0]
    return s(Partition(lam)/Partition(mu))

#Sc(snake matroid)
def ScSnake(b, func = True, eta = None, comp = False):
    k = len(b)
    n = sum(b) + 1
    symfunc = ribbonSchur(b) # = Sc^c(S(b))
    if eta != None:
        if comp: return symfunc.coefficient(eta)
        else: return symfunc.coefficient(complement(k,n,eta))
    if func:
        if comp:
            return symfunc
        else:
            return sum([coeff*s(complement(k,n,eta)) for (eta, coeff) in symfunc])
    if comp: return {tuple(etac) : symfunc.coefficient(etac) for etac in partitions(k,n)}
    else: return {tuple(complement(k,n,etac)) : symfunc.coefficient(etac) for etac in partitions(k,n)}

#list of compositions corresponding to ribbons in lam/mu.
def snakeinLPM(lam, mu):
    k = len(lam)
    n = k + lam[0]
    mu += [0]*(k-len(mu))
    right = [lam[k-i-1]+i for i in range(k-1)]
    left = [mu[k-i-2]+1+i for i in range(k-1)]
    valid = []
    for b in Compositions(n-1,length = k):
        bsum = [sum(b[:s]) for s in range(1,k)]
        if all([left[i] <= bsum[i] and bsum[i] <= right[i] for i in range(k-1)]):
            valid.append(b)
    return valid

#Sc(LPM(lam/mu)) 
#ex: Sc(U24) = ScLPM([2,2],[0])
def ScLPM(lam, mu, func = True, eta = None, comp = False):
    k = len(lam)
    n = k + lam[0]
    symfunc = sum([ribbonSchur(b) for b in snakeinLPM(lam, mu)])
    if eta != None:
        if comp: return symfunc.coefficient(eta)
        else: return symfunc.coefficient(complement(k,n,eta))
    if func:
        if comp:
            return symfunc
        else:
            return sum([coeff*s(complement(k,n,eta)) for (eta, coeff) in symfunc])
    if comp: return {tuple(etac) : symfunc.coefficient(etac) for etac in partitions(k,n)}
    else: return {tuple(complement(k,n,etac)) : symfunc.coefficient(etac) for etac in partitions(k,n)}

#Sc(N)
#N has chain of cyclic flasts of sizes hi and rank ri. 
#data = [(h1,r1),(h2,r2),...,(hs,rs)]
def ScNested(data, func = True, eta = None, comp = False):
    n = data[-1][0]
    k = data[-1][1]
    lam = [n-k]*k
    mu = [0]
    for i in range(len(data)-1):
        mu += [data[i][0] - data[i][1]]*(data[i+1][1] - data[i][1])
    mu = mu[::-1]
    return ScLPM(lam, mu, func = func, eta = eta, comp = comp)

load('nestedbasis.sage')
def Sc(M, func = True, eta = None, comp = False):
    k = M.rank()
    n = len(M.groundset())
    symfunc = 0 #Sc^c(M)
    for (coeff,data) in nestedBasis(M):
        symfunc += coeff*ScNested(data, func = True, comp = True)
    if eta != None:
        if comp: return symfunc.coefficient(eta)
        else: return symfunc.coefficient(complement(k,n,eta))
    if func:
        if comp:
            return symfunc
        else:
            return sum([coeff*s(complement(k,n,eta)) for (eta, coeff) in symfunc])
    if comp: return {tuple(etac) : symfunc.coefficient(etac) for etac in partitions(k,n)}
    else: return {tuple(complement(k,n,etac)) : symfunc.coefficient(etac) for etac in partitions(k,n)}


def example():
    S = snakeMatroid((2,1,2,3))
    print('All possible outputs for S(2,1,2,3)')
    print(ScSnake((2,1,2,3)), '\n')
    print(ScSnake((2,1,2,3), comp = True), '\n')
    print(ScSnake((2,1,2,3), func = False), '\n')
    print(ScSnake((2,1,2,3), func = False, comp = True), '\n')
    print('d_[4,4,3,1] = ', ScSnake((2,1,2,3), eta = [4,4,3,1]), '\n')
    print('d_[4,2,1,1]^c = ', ScSnake((2,1,2,3), eta = [4,2,1,1], comp = True), '\n')

    #snake and nested matroid as LPM.
    print(ScSnake((2,1,2,3)) == ScLPM([5,3,2,2],[2,1,1]))
    print(ScNested([(2,1),(4,2),(7,3)]) == ScLPM([4,4,4],[2,1]), '\n')

    F = matroids.catalog.Fano()
    print('Sc(Fano) = ', Sc(F))



    
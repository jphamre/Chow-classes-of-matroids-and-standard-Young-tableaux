#Functions on matroids:

Ring = PolynomialRing(ZZ,'t')
t = Ring.0

#characteristic polynomial
def charpoly(M):
    return (-1)^(M.rank()) * M.tutte_polynomial(1-t,0)+t-t

#beta invariant
def b(M):
    return (-1)^(M.rank()-1)*diff(charpoly(M),t)(1)

#Operations on matroids:

#Change indexation of matroid by j.
def reindex(M, j = 1):
    EE = [i + j for i in M.groundset()]
    BB = [[i+j for i in B] for B in M.bases()]
    return Matroid(groundset = EE, bases = BB)

#Change indexation of matroid to start at 1.
def oneindex(M):
    j = 1 - min(M.groundset())
    return reindex(M, j)

#Ordered direct sum of matroids. 
def directsum(M1,M2):
    if len(M1.groundset()) == 0:
        return M2
    if len(M2.groundset()) == 0:
        return M1
    M1 = oneindex(M1)
    M2 = oneindex(M2)
    E1 = M1.groundset()
    n = len(E1)
    M2 = reindex(M2, n)
    E2 = M2.groundset()
    EE = E1.union(E2)
    BB = [p[0].union(p[1]) for p in cartesian_product([list(M1.bases()), list(M2.bases())])]
    return Matroid(groundset = EE, bases = BB)

#Extend M with element in parallel with n
def parext(M):
    n = list(M.groundset())[-1]
    return M.extension(n+1,[frozenset([n])])

#Extend M with element in series with n
def serext(M):
    return parext(M.dual()).dual()

#Restriction of M to subset F of [n].
def restrict(M, F):
    E = M.groundset()
    C = E.difference(F)
    return M.delete(C)

#Constructions of matroids:

#Transversal matroids for set system A
def transversalMatroid(setsystem):
    CP = [Set(list(s)) for s in cartesian_product(setsystem)]
    return Matroid(bases = [s for s in CP if len(s) == len(setsystem)])

#LPM by indices of north steps of upper and lower path.
def LPMpaths(U,L):
    setsystem = [list(range(U[i], L[i] + 1)) for i in range(len(L))]
    return transversalMatroid(setsystem)

#Indeces of north steps of lower path of lam.
def lowerpath(lam):
    k = len(lam)
    return [lam[k-1-i] + i + 1 for i in range(k)]

#Connected LPM by skew tableau lam/mu.
def LPM(lam, mu):
    L = lowerpath(lam)
    mu = [i for i in mu if i != 0]
    dif = len(lam) - len(mu)
    U = list(range(1,dif+1)) + [n + dif for n in lowerpath(mu)]
    return LPMpaths(U,L)

#Snake matroid associated to composition b.
def snakeMatroid(b):
    L = [sum(b[:i]) + 1 for i in range(1,len(b)+1)]
    U = [1] + L[:-1]
    return(LPMpaths(U,L))

#A nested matroid with cyclic flats of size hi and rank ri.
# data = [(h1,r1),(h2,r2),...,(n,r)]
def nestedMatroid(data):
    n = data[-1][0]
    k = data[-1][1]
    lam = [n-k]*k
    mu = [0]
    for i in range(len(data)-1):
        mu += [data[i][0] - data[i][1]]*(data[i+1][1] - data[i][1])
    mu = mu[::-1]
    return LPM(lam,mu)

#minimal matroid T24
#= transversalMatroid([[1,2,3],[3,4]])
#= LPMpaths([1,3],[3,4])
#= LPM([2,2],[1])
#= snakeMatroid([2,1])
#= nestedMatroid([(2,1),(4,2)])
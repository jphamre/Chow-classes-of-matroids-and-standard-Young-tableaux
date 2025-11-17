
#Following 'Intersection ring of matroids' - Simon Hampe.
load('operations.sage')

#Lattice of cyclic flats, with rank data.
def lcf(M):
    CF = []
    for F in M.lattice_of_flats():
        if(len(restrict(M,F).coloops()) == 0):
            CF.append(F)
    P = Poset((CF,lambda p,q : p.issubset(q)))
    return P.relabel(lambda cf : (cf, M.rank(cf)))

#Lattice of chains in P with artificial one.
def chainposet(P):
    chains = []
    for C in P.chains():
        if (P.top() in C) and (P.bottom() in C):
            chains.append(tuple(C))
    chains.append('one')
    def order(p,q):
        if p == 'one':
            return q == 'one'
        if q == 'one':
            return True
        return frozenset(p).issubset(frozenset(q))
    return Poset((chains, order))

def isoclass(data):
    return [(len(data[i][0]), data[i][1]) for i in range(1, len(data))]

#Nested basis up to isomorphism.
def nestedBasis(M, matroid = False):
    nested = {} #nested[N] = sign
    cp = chainposet(lcf(M))
    for data in cp:
        if data != 'one':
            coeff = -cp.moebius_function(data ,'one')
            data = tuple(isoclass(data))
            if data in nested:
                nested[data] += coeff
            else:
                nested[data] = coeff
    if matroid:
        return [(nested[data], nestedMatroid(data)) for data in nested.keys()]
    return [(nested[data], data) for data in nested.keys()]
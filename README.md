# Chow-classes-of-matroids-and-standard-Young-tableaux
Compute the class of a matroid in the Chow ring of the Grassmannian as in https://arxiv.org/abs/2511.01711.

## operations.sage
Contains common operations on matroids, and functions to create different types of matroids with our convensions.

## nestedbasis.sage
Computes the coefficients of a matroid in the basis for valuativity of nested matroids, see https://www.sciencedirect.com/science/article/pii/S0095895616300648.

## chowclass.sage
Computing the chow class of matroids. There are spesific functions for different classes of matroids:
- snake matroids, ScSnake
- connected lattice path matroids, ScLPM
- nested matroids, ScNested
- any matroid using nestedbasis, Sc.

Each function has the same options for different types of outputs (see the example function):
- Sc(M) = Sc(M) as a symmetric function
- Sc^c(M) = Sc(M, comp = True) as a symmetric function
- Sc(M) = Sc(M, func = False) as a dictionary
- Sc^c(M) = Sc(M, func = False, comp = True) as a dictionary
- d_eta(M) = Sc(M, eta = [eta1,eta2,...]) an integer
- d_{eta^c}(M) = Sc(M eta = [eta1,eta2,...], comp = True) an integer
  

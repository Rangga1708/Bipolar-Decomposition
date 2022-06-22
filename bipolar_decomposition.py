import numpy as np
import scipy.linalg as la

def diagonalize_orthogonally(A):
    eigenvalues, eigenvectors = la.eig(A)
    
    #orthonormalization (Gram-Schmidt)
    for i in range(A.shape[0]):
        eigenvectors[:,i]=(1/(np.sqrt(np.dot(eigenvectors[:,i],np.conj(eigenvectors[:,i]))))) * eigenvectors[:,i]
    
    return(np.matrix(np.diag(eigenvalues)), np.matrix(eigenvectors))

def bipolar(Z):
    # check if Z is nonsingular matrix
    if abs(la.det(Z)) < 1e-5:
        raise NameError("The matrix is singular. Can't use bipolar decomposition.")

    # calculate A#conj(A), where A = Z*Z
    # find S such that A#conj(A) = e^{2S}
    # we can easily get S = 1/2 log(A#conj(A))
    A = Z.getH() * Z
    A_Aconj = A * la.sqrtm(la.inv(A) * np.conj(A))
    S = np.matrix(np.real(1/2 * la.logm(A_Aconj)))
    
    # calculate Y = e^{-S}Ae^{-S} dan find K such that Y = e^{2iK}
    # we can easily get iK = 1/2 log(Y)
    Y = la.expm(-S) * A * la.expm(-S)
    K = np.matrix(np.imag(1/2 * la.logm(Y)))

    # calculate W
    # calculate W'W dan find T such that W'W = e^{2iT}
    # we can easily get iT = 1/2 log(W'W)
    W = Z * la.expm(-S) * la.expm(-1j * K)
    Wtrans_W = W.getT() * W
    T = np.matrix(np.imag(1/2 * la.logm(Wtrans_W)))
    
    # calculate Q = We^{-iT}
    # if det(Q) = 1, find L such that Q = e^L
    # if det(Q) = -1, we need further steps to find L
    Q = np.matrix(np.real(W * la.expm(-1j * T)))
    if la.det(Q) > 0:
        L = np.matrix(la.logm(Q))
    else:
        # first, diagonalize T
        _, G = diagonalize_orthogonally(T)
        
        # calculate X = GDG', where D = diag(pi,0,...,0)
        # calculate J = e^{iX}
        D = np.zeros(Z.shape)
        D[0,0] = np.pi
        X = G * D * G.getT()
        J = np.real(la.expm(1j * X))
        
        # calculate QJ and find L such that QJ = e^L
        # relabel T = X + T
        QJ = np.real(Q * J)
        L = np.matrix(la.logm(QJ))
        T = X + T
        
    return(L,T,K,S)
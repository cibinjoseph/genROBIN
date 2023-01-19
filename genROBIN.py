import numpy as np

nxFuselage = 24
ntFuselage = 16
nxPylon = 12
ntPylon = 12

fusFile = "robinFuselage.obj"
pylFile = "robinPylon.obj"


def getSectionIndex(x, isPylon = False):
    """ Returns an index corresponding to each section
        for obtaining values from coeff matrix """
    idx = -1
    if isPylon:
        idx = 4 if x < 0.8 else 5
    else:
        if x < 0.4:
            idx = 0
        elif x < 0.8:
            idx = 1
        elif x < 1.9:
            idx = 2
        else:
            idx = 3
    return idx

def getChebyshevNode(a, b, k, n):
    if k > 0 and k < n:
        return 0.5*(a+b) + 0.5*(b-a)*np.cos((2.0*(n-k))*np.pi*0.5/n)
    elif k == 0:
        return a
    else:
        return b

def createVertices(nx, nt, filename, isPylon = False):
    """ Creates vertices for fuselage or pylon geometry """

    # Rows 0, 1, 2 and 3 of the coefficient matrix are for the fuselage
    # Rows 4 and 5 are for the pylon
    hcoeff = np.array([[1.0, -1.0, -0.4, -0.4,   1.8, 0.0,  0.25,  1.8],
                       [0.0,  0.0,  0.0,  1.0,   0.0, 0.25, 0.0,   1.0],
                       [1.0, -1.0, -0.8,  1.1,   1.5, 0.05, 0.2,   0.6],
                       [1.0, -1.0, -1.9,  0.1,   2.0, 0.0,  0.05,  2.0],
                       [1.0, -1.0, -0.8, -0.4,   3.0, 0.0,  0.145, 3.0],
                       [1.0, -1.0, -0.8,  0.218, 2.0, 0.0,  0.145, 2.0]])

    wcoeff = np.array([[1.0, -1.0, -0.4, -0.4,   2.0, 0.0,  0.25,  2.0],
                       [0.0,  0.0,  0.0,  1.0,   0.0, 0.25, 0.0,   1.0],
                       [1.0, -1.0, -0.8,  1.1,   1.5, 0.05, 0.2,   0.6],
                       [1.0, -1.0, -1.9,  0.1,   2.0, 0.0,  0.05,  2.0],
                       [1.0, -1.0, -0.8, -0.4,   3.0, 0.0,  0.166, 3.0],
                       [1.0, -1.0, -0.8,  0.218, 2.0, 0.0,  0.166, 2.0]])

    zcoeff = np.array([[1.0, -1.0, -0.4, -0.4, 1.8, -0.08,  0.08, 1.8],
                       [0.0,  0.0,  0.0,  1.0, 0.0,  0.0,   0.0,  1.0],
                       [1.0, -1.0, -0.8,  1.1, 1.5,  0.04, -0.04, 0.6],
                       [0.0,  0.0,  0.0,  1.0, 0.0,  0.04,  0.0,  1.0],
                       [0.0,  0.0,  0.0,  1.0, 0.0,  0.125, 0.0,  1.0],
                       [1.0, -1.0, -0.8,  1.1, 1.5,  0.065, 0.06, 0.6]])

    ncoeff = np.array([[2.0,  3.0,  0.0, 0.4, 1.0, 0.0, 1.0, 1.0],
                       [0.0,  0.0,  0.0, 1.0, 0.0, 5.0, 0.0, 1.0],
                       [5.0, -3.0, -0.8, 1.1, 1.0, 0.0, 1.0, 1.0],
                       [0.0,  0.0,  0.0, 1.0, 0.0, 2.0, 0.0, 1.0],
                       [0.0,  0.0,  0.0, 1.0, 0.0, 5.0, 0.0, 1.0],
                       [0.0,  0.0,  0.0, 1.0, 0.0, 5.0, 0.0, 1.0]])

    # Fixes from Applied-Scientific-Research/robin-surface-mesh.git
    # 1) if there's a 0.0 in the second col, then change the 4th and 5th cols to 1.0
    # 2) if there's a 0.0 in C7, change C8 to 1.0, same as above, to prevent nan/inf
    # 3) the 0.4..0.8 section (row 2) coefficients in C1 needed to go into C6
    # 4) C4 is wrong in the first section of fuse and pyl - it needed to be negative

    if isPylon:
        # Pylon section limits
        xBegin = 0.4
        xEnd = 1.018
    else:
        # Fuselage section limits
        xBegin = 0.0
        xEnd = 2.0

    with open(filename, "w") as fh:
        fh.write("# Vertices\n")

        for ix in range(nx+1):
            xol = getChebyshevNode(xBegin, xEnd, ix, nx)
            sec = getSectionIndex(xol)

            if sec == -1:
                raise ValueError("Incorrect Chebyshev node value")

    return


if __name__ == "__main__":

    # Check if inputs are valid
    if all((nxFuselage, ntFuselage, nxPylon, ntPylon)) > 0:
        print("Generating ROBIN geometry")
    else:
        raise ValueError("Incorrect number of elements")

    # Create fuselage
    print("Generating fuselage")
    createVertices(nxFuselage, ntFuselage, fusFile)
    # -- DEBUG --
    exit()
    createFaces(fusFile, nxFuselage, ntFuselage)

    # Create pylon
    print("Generating pylon")
    createVertices(nxFuselage, ntFuselage, fusFile, isPylon = True)
    createFaces(fusFile, nxFuselage, ntFuselage, isPylon = True)

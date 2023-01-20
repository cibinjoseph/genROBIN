import numpy as np

nxFuselage = 12
ntFuselage = 8
nxPylon = 6
ntPylon = 6

fusFile = "robinFuselage.obj"
pylFile = "robinPylon.obj"

eps = np.spacing(1.0)

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

def getsuperval(x, c):
    val = c[5] + c[6]*pow(max(0.0, c[0] + \
                 c[1]*pow((x+c[2])/c[3], c[4])), 1.0/c[7])
    return val

def getRadialCoordinate(H, W, theta, N):
    numer = 0.25*H*W
    denom = pow(0.5*H*np.abs(np.sin(theta)), N) + \
            pow(0.5*W*np.abs(np.cos(theta)), N)
    if (abs(denom) < eps):
        denom = 1
    return numer / pow(denom, 1.0/N)

def getVertices(nx, nt, isPylon = False):
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

    xol = np.empty(shape=[nx+1, nt])
    yol = np.empty(shape=[nx+1, nt])
    zol = np.empty(shape=[nx+1, nt])

    for ix in range(nx+1):
        xval = getChebyshevNode(xBegin, xEnd, ix, nx)
        xol[ix, :] = xval
        sec = getSectionIndex(xval)

        if sec == -1:
            raise ValueError("Incorrect Chebyshev node value")

        H  = getsuperval(xval, hcoeff[sec, :])
        W  = getsuperval(xval, wcoeff[sec, :])
        Z0 = getsuperval(xval, zcoeff[sec, :])
        N  = getsuperval(xval, ncoeff[sec, :])

        for it in range(nt):
            theta = 2.0*np.pi*it/float(nt)

            r = getRadialCoordinate(H, W, theta, N)
            yol[ix, it] = r * np.sin(theta)
            zol[ix, it] = r * np.cos(theta) + Z0

            if (ix == 0) or (ix == nx):
                # Avoid multiple coincident points at
                # tip and tail
                break
    return xol, yol, zol

def writeVertices(xol, yol, zol, filename):
    nx, nt = xol.shape
    with open(filename, "w") as fh:
        fh.write("# Vertices\n")
        for ix in range(nx):
            for it in range(nt):
                fh.write("v " + \
                         str(xol[ix, it]) + " " + \
                         str(yol[ix, it]) + " " + \
                         str(zol[ix, it]) + "\n")
                if (ix == 0) or (ix == nx-1):
                    # Avoid multiple coincident points at
                    # tip and tail
                    break

if __name__ == "__main__":

    # Check if inputs are valid
    if all((nxFuselage, ntFuselage, nxPylon, ntPylon)) > 0:
        print("Generating ROBIN geometry")
    else:
        raise ValueError("Incorrect number of elements")

    # Create fuselage
    print("Generating fuselage")
    x, y, z = getVertices(nxFuselage, ntFuselage)
    writeVertices(x, y, z, fusFile)
    # DEBUG
    exit()
    createFaces(fusFile, nxFuselage, ntFuselage)

    # Create pylon
    print("Generating pylon")
    x, y, z = createVertices(nxFuselage, ntFuselage, isPylon = True)
    createFaces(fusFile, nxFuselage, ntFuselage, isPylon = True)

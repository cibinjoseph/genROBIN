#!/usr/bin/python3

import numpy as np
import argparse as ag

# Machine epsilon
eps = np.spacing(1.0)

def getSectionIndex(x, isPylon = False):
    """ Returns an index corresponding to each section
        for obtaining values from coeff matrix """
    idx = np.empty(x.size).astype(int)
    for i in range(x.size):
        if isPylon:
            idx[i] = 4 if x[i] < 0.8 else 5
        else:
            if x[i] < 0.4:
                idx[i] = 0
            elif x[i] < 0.8:
                idx[i] = 1
            elif x[i] < 1.9:
                idx[i] = 2
            else:
                idx[i] = 3
    return idx

def getChebyshevNodes(a, b, n):
    """ Get n+1 Chebyshev nodes from a to b """
    k = np.arange(n+1)
    nodes = 0.5*(a+b) + 0.5*(b-a)*np.cos((2.0*(n-k))*np.pi*0.5/n)
    nodes[0] = a
    nodes[n] = b
    return nodes

def getsuperval(x, c):
    cval = (x+c[2])/c[3]
    # np.power(x) does not like negative terms
    # Computing using workaround np.sign(x)*np.abs(x)**y
    negPowerTerm = c[0] +c[1]*np.sign(cval)*np.abs(cval)**c[4]
    val = c[5] + c[6]*np.power(np.maximum(0.0, negPowerTerm), 1.0/c[7])
    return val

def getRadialCoordinate(H, W, theta, N):
    """ Returns radial coordinate for the np.array theta """
    numer = 0.25*H*W
    denom = np.power(0.5*H*np.abs(np.sin(theta)), N) + \
            np.power(0.5*W*np.abs(np.cos(theta)), N)
    # Set zero elements in denominator to 1
    denom[np.abs(denom) < eps] = 1.0
    return numer / np.power(denom, 1.0/N)

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

    xval = getChebyshevNodes(xBegin, xEnd, nx)
    xol = np.tile(xval, (nt, 1)).T

    secIdx = getSectionIndex(xval, isPylon)
    theta = 2.0*np.pi*np.arange(nt)/float(nt)

    for ix in range(nx+1):
        H  = getsuperval(xval, hcoeff[secIdx[ix], :])
        W  = getsuperval(xval, wcoeff[secIdx[ix], :])
        Z0 = getsuperval(xval, zcoeff[secIdx[ix], :])
        N  = getsuperval(xval, ncoeff[secIdx[ix], :])

        r = getRadialCoordinate(H[ix], W[ix], theta, N[ix])
        yol[ix, :] = np.multiply(r, np.sin(theta))
        zol[ix, :] = np.multiply(r, np.cos(theta)) + Z0[ix]

    return xol, yol, zol

def getFaceInfo(x, y, z):
    """ Creates face information for triangle surface elements """
    nx1, nt = x.shape
    nx = nx1 - 1

    sz = 2*nt*(nx-1)
    faceNodes = np.empty(shape=[sz, 3]).astype(int)

    idx = 0
    for i in range(2, nt+1):
        faceNodes[idx, :] = [1, i, i+1]
        idx += 1
    faceNodes[idx, :] = [1, nt+1, 2]
    idx += 1

    for r in range(nx-2):
        ir1 = nt*r + 2
        ir2 = nt*(r+1) + 2
        for n in range(nt-1):
            faceNodes[idx, :] = [ir1+n, ir2+n, ir1+n+1]
            idx += 1
            faceNodes[idx, :] = [ir2+n, ir2+n+1, ir1+n+1]
            idx += 1
        faceNodes[idx, :] = [ir1+nt-1, ir2+nt-1, ir1]
        idx += 1
        faceNodes[idx, :] = [ir2+nt-1, ir2, ir1]
        idx += 1

    lastVert = nt*(nx-1)+2
    for i in range(2, nt+1):
        faceNodes[idx, :] = [(i+1)+nt*(nx-2) ,i+nt*(nx-2), lastVert]
        idx += 1

    faceNodes[idx, :] = [2+nt*(nx-2), lastVert-1, lastVert]

    return faceNodes

def writeVertices(xol, yol, zol, filename, writeFaces=True):
    """ Writes vertices and faces to file """
    nx1, nt = xol.shape
    nx = nx1 - 1

    with open(filename, "w") as fh:
        fh.write("# Vertices\n")
        for ix in range(nx+1):
            for it in range(nt):
                fh.write(f"v {xol[ix, it]} {yol[ix, it]} {zol[ix, it]}\n")
                if (ix == 0) or (ix == nx):
                    # Avoid multiple coincident points at
                    # tip and tail
                    break

        if writeFaces:
            fh.write("\n# Faces\n")
            face = getFaceInfo(x, y, z)
            for i in range(face.shape[0]):
                fh.write(f"f {face[i, 0]} {face[i, 1]} {face[i, 2]}\n")

def getArguments():
    # Create parser for arguments
    parser = ag.ArgumentParser()
    parser.add_argument("nxFuselage", type=int, \
                        help="No. of lengthwise elements for fuselage")
    parser.add_argument("ntFuselage", type=int, \
                        help="No. of circumferential elements for fuselage")
    parser.add_argument("nxPylon", type=int, \
                        help="No. of lengthwise elements for pylons")
    parser.add_argument("ntPylon", type=int, \
                        help="No. of circumferential elements for pylons")

    return parser.parse_args()

if __name__ == "__main__":

    # Parse arguments
    args = getArguments()

    nscale = 1
    nxFuselage = nscale*args.nxFuselage
    ntFuselage = nscale*args.ntFuselage
    nxPylon = nscale*args.nxPylon
    ntPylon = nscale*args.ntPylon

    fileType = "obj"

    fusFile = "robinFuselage." + fileType
    pylFile = "robinPylon." + fileType

    # Check if inputs are valid
    if all((nxFuselage, ntFuselage, nxPylon, ntPylon)) > 0:
        print("Generating ROBIN geometry")
    else:
        raise ValueError("Incorrect number of elements")

    # Create fuselage
    x, y, z = getVertices(nxFuselage, ntFuselage)
    print("Writing fuselage to " + fusFile)
    writeVertices(x, y, z, fusFile, writeFaces=True)

    # Create pylon
    x, y, z = getVertices(nxPylon, ntPylon, isPylon = True)
    print("Writing pylon to " + pylFile)
    writeVertices(x, y, z, pylFile, writeFaces=True)

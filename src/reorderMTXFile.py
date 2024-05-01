import sys
import numpy as np
import scipy


""" reorderMTXFile.py - Script for converting MTX format
matrices from column-wise ordering to row-wise ordering.

This python code will replace col-wise ordering to row-wise ordering.

Usage:
    python3 reorderMTXFile.py <input_filename>

Example:
    python3 reorderMTXFile.py cage3.mtx

This will write the reordered matrix data to a new file named cage3_T.mtx.
"""



def readMTXFile(filename):
    F = open(filename, "r")
    line = F.readline()

    preLines = list()
    # Step 1: Ignoring first lines. Writing
    while line:
        if line[0] != '%':
            break
        preLines.append(line)
        line = F.readline()


    # Step 2: Get number of rows, col, and number of non-zero elements.
    N, M, NNZ = line.split()
    N, M, NNZ = int(N), int(M), int(NNZ)
    print(f"N, M, NNZ = {N}, {M}, {NNZ}")

    row = np.zeros(NNZ, dtype=np.int32)
    col = np.zeros(NNZ, dtype=np.int32)
    val = np.zeros(NNZ)


    # Step 3: Get the values of the matrix.
    line = F.readline()
    nnz = 0;
    while line:
        i, j, x = line.split()
        i, j, x = int(i)-1, int(j)-1, float(x)

        row[nnz] = i
        col[nnz] = j
        val[nnz] = x

        nnz += 1
        line = F.readline()
    F.close()

    A = scipy.sparse.csr_matrix((val, (row, col)))
    # print(A.toarray())
    return A, NNZ, preLines


def writeMatrixToMTX(A, NNZ, filename, preLines):
    """ Writes a matrix to a file in MTX format.
    Using row-wise ordering instead of col-wise. """

    F = open(filename, "w")
    for line in preLines:
        F.write(line)
    N, M = A.shape
    print(f"N, M = {N}, {M}")

    line = str(N) + " " + str(M) + " " + str(NNZ) + "\n"
    F.write(line)

    rows, cols = A.nonzero()
    for row,col in zip(rows, cols):
        # print(f"(row, col), A[row,col] = ({row}, {col}), {A[row, col]}")
        line = str(row+1) + " " + str(col+1) + " " + str(A[row, col]) + "\n"
        F.write(line)


    F.close()
    return


def main():
    filename_pre = "/mnt/c/Users/xu/petsc/share/petsc/datafiles/matrices/MYMAT/"
    filename_mid = "1138_bus"
    if len(sys.argv) >= 2:
        filename_mid = sys.argv[1]
    filename_post = ".mtx"

    inputFile = filename_pre + filename_mid + filename_post
    outputFile = filename_pre + filename_mid + "_T" + filename_post
    print(f"input filename = {inputFile}")
    print(f"outputfilename = {outputFile}\n")


    A, NNZ, context = readMTXFile(inputFile)

    if (A.getnnz() < 25):
        print(A.toarray())
        print()

    writeMatrixToMTX(A, NNZ, outputFile, context)











if __name__ == "__main__":
    main()

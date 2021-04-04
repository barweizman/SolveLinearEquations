# ------Solve AX=b by using Gauss Elimination------

# Singularity check-by calculating determination
def CalcDet(matrix):
    """
    description:
    recursive function.calculate the determination of the matrix
    :param matrix: array nXn. represent the matrix
    :return:matrix's determination
    """
    n = len(matrix)
    if n == 1:
        return matrix[0][0]
    det = 0
    for j in range(0, n):
        det += ((-1) ** j) * matrix[0][j] * CalcDet(SubDet(matrix, 0, j))
    return det


def SubDet(matrix, i, j):
    """
    description: the function create new matrix (n-1)(n-1) by
    removing the selected column and row.
    :param i: row to remove
    :param matrix: the original matrix
    :param j:column to remove
    :return: minor matrix
    """
    return [row[:j] + row[j + 1:] for row in (matrix[:i] + matrix[i + 1:])]


# Elementary Actions on Matrix-Creating an elementary matrix
def Identity(n):
    """
    description:
    The function Create identity matrix
    :param n: size of choice
    :return: identity matrix
    """
    mat = [([0] * n) for i in range(n)]  # initialize the matrix with zeros
    for i in range(0, n):
        mat[i][i] = 1  # the identity matrix includes 1 all over its diagonal, starts at [0][0]
    return mat


def Matrix_multiplication(mat1, mat2):
    """
    description:
    The function the result matrix of multiplication two matrix
    mat1 on the left and mat2 on the right
    :param mat1: first initialized square matrix with values
    :param mat2: second initialized square matrix with values
    :return: result matrix, that will be the multiplication between mat1 and mat2
    """
    if len(mat1[0]) != len(mat2):
        raise Exception("Illegal multiplication between matrix's ")
    result_mat = [([0] * len(mat2[0])) for i in range(len(mat1))]  # initialize the result matrix with zeros

    # iterate through the first matrix rows
    for row1 in range(0, len(mat1)):
        # iterate through the second matrix columns
        for col2 in range(0, len(mat2[0])):
            # iterate through the second matrix rows
            for row2 in range(0, len(mat2)):
                result_mat[row1][col2] += mat1[row1][row2] * mat2[row2][col2]
    return result_mat


def Matrix_addition(mat1, mat2):
    """
    :param mat1: first initialized square matrix with values
    :param mat2: second initialized square matrix with values
    :return: result matrix, that will be the addition between mat1 and mat2
    """
    size = len(
        mat1)  # doesnt matter what we choose(mat1,mat2,rows or columns), we take into consideration only square matrix
    result_mat = [([0] * size) for i in range(size)]  # initialize the result matrix with zeros

    # iterate through the first matrix rows
    for row in range(0, len(mat1)):
        # iterate through the first matrix columns
        for col in range(0, len(mat1[0])):
            result_mat[row][col] = mat1[row][col] + mat2[row][col]
    return result_mat


def Copy_matrix(mat):
    """
    description:
    the function receive matrix and return a copy of her
    :param mat: matrix with values
    :return: copy of the given matrix
    """
    size = len(mat)
    copied_mat = [([0] * size) for i in range(size)]  # initialize the matrix with zeros
    for row in range(0, len(mat)):
        for col in range(0, len(mat[0])):
            copied_mat[row][col] = mat[row][col]  # copy all values
    return copied_mat


def ResetOrgan(row, col, n, pivot, a):
    """
    description:
    create elementary matrix that reset the chosen organ to zero
    at first the elementary matrix starts as a regular I matrix.
    at the chosen indexes we put the value which will give a 0 by multiply the elementary matrix with the original matrix
    :param row: num of row to insert the value
    :param col: num of column to insert the value
    :param n: number of rows and column
    :param pivot: the pivot of the original matrix
    :param a: the value we need to reset to zero
    :return:elementary matrix
    """
    elementary_matrix = Identity(n)
    elementary_matrix[row][col] = -(a / pivot)
    return elementary_matrix


def MultiplyRow(row, a, n):
    """
    description:
    The function create elementary matrix through which we can multiply a row by a value
    insert the value to the pivot in the selected row
    :param row:the selected row to multiply
    :param a:the value to multiply the row
    :param n:number of rows and columns in the original matrix
    :return:elementary matrix
    """
    elementary_matrix = Identity(n)
    elementary_matrix[row][row] = a
    return elementary_matrix


def ExchangeRows(row1, row2, n):
    """
    The function creates for us an elementary matrix through which we can exchange rows
    by multiplying them
    :param row1:row to exchange
    :param row2:row to exchange
    :param n:number of rows and columns in the original matrix
    :return:elementary matrix
    """
    elementary_matrix = Identity(n)
    elementary_matrix[row1][row1] = 0
    elementary_matrix[row1][row2] = 1
    elementary_matrix[row2][row2] = 0
    elementary_matrix[row2][row1] = 1
    return elementary_matrix


def InvertMatrix(matrix):
    """
    description:
    The function calculate the inverted matrix. suitable for singular matrix-if not singular an exception will be raised
    :param matrix: the matrix to invert
    :return: inverted matrix
    """
    if len(matrix) != len(matrix[0]):
        raise Exception("not a singular matrix. there is no inverted matrix")
    n = len(matrix)
    result = Identity(n)
    for j in range(0, n):
        for i in range(0, n):
            if i == j:
                if matrix[i][j] != 0:
                    pivot = matrix[i][j]
                else:
                    for k in range(i + 1, n):
                        if matrix[k][j] != 0:
                            elementarymatrix = ExchangeRows(k, i, n)
                            matrix = Matrix_multiplication(elementarymatrix, matrix)
                            result = Matrix_multiplication(elementarymatrix, result)
                            pivot = matrix[i][j]
                            break
                    if matrix[i][j] == 0:
                        raise Exception("not a singular matrix. there is no inverted matrix")

        for i in range(0, n):
            if i != j:
                if matrix[i][j] != 0:
                    elementarymatrix = ResetOrgan(i, j, n, pivot, matrix[i][j])
                    matrix = Matrix_multiplication(elementarymatrix, matrix)
                    result = Matrix_multiplication(elementarymatrix, result)

    for i in range(0, n):
        if matrix[i][i] != 1:
            if matrix[i][i] < 0:
                elementarymatrix = MultiplyRow(i, -1, n)
                matrix = Matrix_multiplication(elementarymatrix, matrix)
                result = Matrix_multiplication(elementarymatrix, result)

            elementarymatrix = MultiplyRow(i, 1 / matrix[i][i], n)
            matrix = Matrix_multiplication(elementarymatrix, matrix)
            result = Matrix_multiplication(elementarymatrix, result)
    for row in range(n):
        for col in range(n):
            result[row][col] = round(result[row][col], 2)
    return result


# two ways to calculate matrix
def LU(mat):
    """

    :param mat:
    :return:
    """
    size = len(mat)  # doesnt matter if we take len(mat) or len(mat[0]), were talking square matrix.
    res = Copy_matrix(mat)
    L = [([0] * size) for i in range(size)]  # initialize the matrix with zeros
    for col in range(0, size):
        pivot = mat[col][col]
        for row in range(0, size):
            m = -(mat[row][col] / pivot)  # m is the multiply number
            elementary_mat = Identity(size)
            elementary_mat[row][col] = m  # this is the elementary matrix (the identity matrix and the multiply number)
            res = Matrix_multiplication(elementary_mat, res)  # at the end of the loops, res will be 'U'
            L = Matrix_addition(InvertMatrix(res), L)

    if Matrix_multiplication(L, res) == mat:
        pass


def CalcSingularMatrix(matrix, b):
    """
    description:
    the function calculate the result by multiply the inverted matrix and the result vector
    :param matrix: the given matrix
    :param b: the result vector
    :return:the Vector variables
    """
    return Matrix_multiplication(InvertMatrix(matrix), b)


# main function to calculate matrix
def CalcMatrix(matrix, b):
    if CalcDet(matrix) == 0:
        CalcSingularMatrix(matrix, b)
    else:
        pass


# driver
mat1 = [[2, -1, 0],
        [-1, 2, -1],
        [0, -1, 2]]
mat2 = [[4], [6], [8]]
try:
    # print(InvertMatrix(mat1))
    print(CalcSingularMatrix(mat1, mat2))
except Exception as e:
    print(e)

mat1 = [[1, 0],
        [-1, 3]]

mat2 = [[3, 1],
        [2, 1]]

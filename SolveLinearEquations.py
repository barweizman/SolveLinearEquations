
mat1 = [[1,0],
       [-1,3]]

mat2 = [[3,1],
       [2,1]]


def Matrix_multiplication(mat1, mat2):
    """
    :param mat1: first initialized square matrix with values
    :param mat2: second initialized square matrix with values
    :return: result matrix, that will be the multiplication between mat1 and mat2
    """
    size = len(mat1)  # doesnt matter what we choose(mat1,mat2,rows or columns), we take into consideration only square matrix
    result_mat = [([0] * size) for i in range(size)]   # initialize the result matrix with zeros

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
    size = len(mat1)  # doesnt matter what we choose(mat1,mat2,rows or columns), we take into consideration only square matrix
    result_mat = [([0] * size) for i in range(size)]   # initialize the result matrix with zeros

    # iterate through the first matrix rows
    for row in range(0, len(mat1)):
        # iterate through the first matrix columns
        for col in range(0, len(mat1[0])):
            result_mat[row][col] = mat1[row][col] + mat2[row][col]
    return result_mat

def Identity(n):
    """
    :param n: size of choice
    :return: identity matrix
    """
    mat=[([0] * n) for i in range(n)]  # initialize the matrix with zeros
    for i in range(0,n):
        mat[i][i] = 1  # the identity matrix includes 1 all over its diagonal, starts at [0][0]
    return mat

def Copy_matrix(mat):
    """
    :param mat: matrix with values
    :return: copy of the given matrix
    """
    size = len(mat)
    copied_mat = [([0] * size) for i in range(size)]  # initialize the matrix with zeros
    for row in range(0,len(mat)):
        for col in range(0,len(mat[0])):
            copied_mat[row][col] = mat[row][col]  # copy all values
    return copied_mat


def LU(mat):
    size = len(mat)  # doesnt matter if we take len(mat) or len(mat[0]), were talking square matrix.
    res = Copy_matrix(mat)
    L = [([0] * size) for i in range(size)] # initialize the matrix with zeros
    for col in range(0, size):
        pivot = mat[col][col]
        for row in range(0, size):
            m = -(mat[row][col] / pivot)  # m is the multiply number
            elementary_mat = Identity(size)
            elementary_mat[row][col] = m  # this is the elementary matrix (the identity matrix and the multiply number)
            res = Matrix_multiplication(elementary_mat, res)  # at the end of the loops, res will be 'U'
            L = Matrix_addition(Invert(res), L)

    if Matrix_multiplication(L,res) == mat:
        pass




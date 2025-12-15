#
# Program to work in conjunction with flag.cpp - provides inexact rounding.
#
# Carefull! This works  ONLY if all entries in the semidefinite program has integer
# values!!!!
#  - This can be easily achieved by multiplying the objective function
#  - Also don't forget to multiply linear constraints if you use any
#  - If do not multiply the objective, you can provide a scale to this program to  make
#    the constants integers. No need to resolve the CSDP program. But linear constraints
#    are not scaled by this. Only the objective.
#
# This is a sage code for rounding CSDP programs
# It loads the solution, on every block X it computes decomposition
# that looks like X = A^T D A
# then it rounds A and D to Integers by multiplying by 1000000
# So all etries are integers scaled by 1000000^3
# If X is not a square block (slack variables and linear constraints - block
# with entries only on the diagonal), it is left as a row vector (1 x M) matrix.
#
# Moreover, the equations for flag algebras are actually inequalities.
# It makes no sense to round the 'slack variables' so the slack variables
# are all rounded to 0. Of course, the result are then inequalities and
# we just compute the minimum.
#
# Also, some SDP programs are way too big to fit in memory for rounding
# even on BIG machines.  A small workaround is that first, everything is
# rounded and then the rounding is tested and these can be separate  runs.
# Moreover, testing can be done by parts rather than all at once.
#
# Use Howto:
#
#  def round_program(file, scale, sdp_scaled=True):
#
#

"""
Structures:

problem = dict()
problem['num_constraints'] = num_constraints       # total number of constraints - integer
problem['loaded_constraints'] = used_constraints   # if less than all constraints are loaded
problem['num_blocks'] = num_blocks                 # number of blocks on the program
problem['block_sizes'] = block_sizes               # list of block sizes ( list of integers )
problem['blocks'] = blocks                         # list[constriants][blocks]matrix
problem['constants'] = constants                   # list of constants (one per constriant)

"""

import numpy
import pickle  # for saving/loading the rounded matrices

######################################## Global Vagiables ######################

field = Integer
"""Field used in rounded matrices.
We chose arbitrary precision integers. They are fast for
computation. A drawback is that SAGE does NOT like
matrices over Integers as they are not a field. 
We use 2D lists instead.
"""

roundingDenominator = 100000000000000
"""Scale of all rounded variables. 
Since we use integers, we need to carry decimal places
somehow and we do it by thinking of this as a universal
denominator for all.
"""


eigVZero = 0.00000000001
"""Threshold when number is considered zero"""

runParellel = True
"""Some parts can run in parallel if set to True"""

################################ Loading SDP in dat-s form and CSDP solution ###


def read_problem_outline(filename):
    """Loads dat-s file outline.

    Reads basic part of the problem but does not read the equations.
    It is faster than and more memory efficient to load just the
    outline and it is all that is needed when rounding a solution.
    Args:
        filename: String containing path to dat-s file to load
    Returns:
        problem without constraints
    """
    f = open(filename, "r")

    line = f.readline()
    line = line.strip()
    num_constraints = int(line)

    line = f.readline()
    line = line.strip()
    num_blocks = int(line)

    line = f.readline()
    parts = line.split()
    block_sizes = []
    for part in parts:
        block_sizes.append(int(part))

    problem = dict()
    problem["num_constraints"] = num_constraints
    problem["num_blocks"] = num_blocks
    problem["block_sizes"] = block_sizes
    # problem['blocks'] = blocks            # not read
    # problem['constants'] = constants      # not read
    f.close()
    return problem


def rescale_number(x, constants_rescale):
    """Rescaling and testing that rescaling indeed gives an integer"""

    if constants_rescale == 1:
        return field(x)

    x_float = float(x) * constants_rescale
    x_rounded = round(x_float)
    if abs(x_float - x_rounded) > 0.00001:
        print(
            "Rescaling did not work for x=",
            x,
            " result was",
            x_rounded,
            "which is not precise for",
            x_float,
        )
        print("Try a different rescaling of the semidefinite program.")
        sys.exit()

    return field(x_rounded)


# field(int(float(parts[i])*constants_rescale+0.5))


# Reads the whole problem as integer program
# the program has to be integral. This is a hack
# for SDP programs that were solved without scaled objective
# and are expensive to resolve
#
def read_problem_file(
    filename, constants_rescale=1, constraint_from=0, constraint_to=10000000
):
    """Loads dat-s file.

    Reads dat-s fule and performs rescale if needed. Possible to load only some of the
    constraints from the file in order to save memory.
    Args:
        filename: String containing path to dat-s file to load.
        constants_rescale: Rescaling constants as would correspond to rescale the objective
            function at the beginning before running CSDP
        constraint_from: First constraint to be loaded
        constraint_to: Last constraint to be loaded (capped by number of constraints)
    Returns:
        problem structure
    """
    f = open(filename, "r")

    line = f.readline()
    line = line.strip()
    num_constraints = int(line)

    line = f.readline()
    line = line.strip()
    num_blocks = int(line)

    line = f.readline()
    parts = line.split()
    block_sizes = []
    for part in parts:
        block_sizes.append(int(part))

    line = f.readline()
    parts = line.split()

    constants = []
    used_constraints = 0

    for i in range(num_constraints):
        if i >= constraint_from and i < constraint_to:
            constants.append(rescale_number(parts[i], constants_rescale))
            used_constraints += 1

    """
    if constants_rescale == 1:
        for i in xrange(num_constraints):
           if i >= constraint_from and i < constraint_to:
               constants.append(field(parts[i]))
               used_constraints += 1
    else:
        for i in xrange(num_constraints):
           if i >= constraint_from and i < constraint_to:
               constants.append(field(int(float(parts[i])*constants_rescale+0.5)))
               used_constraints += 1
    """
    # blocks that are integer numbers right away
    blocks = []
    for size in block_sizes:
        if size < 0:
            # block = [ matrix(ZZ,1,-size) for x in range(num_constraints+1) ]
            block = [matrix(ZZ, 1, -size) for x in range(used_constraints + 1)]
        else:
            block = [
                matrix(ZZ, size, size, sparse=True) for x in range(used_constraints + 1)
            ]
            # block = [ matrix(ZZ,size,size, sparse=True) for x in range(num_constraints+1) ]
            # block = [ matrix(QQ,size,size) for x in range(num_constraints+1) ]
        blocks.append(block)

    while True:
        line = f.readline()
        if len(line) == 0:
            break
        parts = line.split()
        if len(parts) == 0:
            continue
        constraint = int(parts[0])
        # Read also a bit arround to accommodate for mistakes..
        # Blocks are shifted by one
        if constraint <= constraint_from or constraint > constraint_to:
            continue
        block_num = int(parts[1])
        row = int(parts[2])
        column = int(parts[3])
        value = rescale_number(parts[4], constants_rescale)
        """
        if constants_rescale == 1:
            value = field(parts[4])
        else:
            value = field(float(parts[4])*constants_rescale)
        """
        if block_sizes[block_num - 1] < 0:
            blocks[block_num - 1][constraint - constraint_from][0, column - 1] += value
        else:
            blocks[block_num - 1][constraint - constraint_from][
                row - 1, column - 1
            ] += value
            if row != column:
                blocks[block_num - 1][constraint - constraint_from][
                    column - 1, row - 1
                ] += value

    problem = dict()
    problem["num_constraints"] = num_constraints
    problem["loaded_constraints"] = used_constraints
    problem["num_blocks"] = num_blocks
    problem["block_sizes"] = block_sizes
    problem["blocks"] = blocks
    problem["constants"] = constants
    f.close()
    return problem


# This reads the solution file - but only the primal. Dual is being ignored
# and it would waste space anyway. The problem is read as real numbers - because
# we
def read_solution_file(filename, problem):
    """Reads solution file (no rounding).
    It reads only the primal from the solution as that is the part for rounding
    and it saves some memory to avoid the dual. The solution is read
    as RDF (real float) numbers.
    Args:
        filename: Name of file containing the solution to be loaded
        problem: problem OUTLINE that corresponds to the solution
    Returns:
        Solution to the primal in RDF matrix form
    """
    f = open(filename, "r")

    num_constraints = problem["num_constraints"]
    num_blocks = problem["num_blocks"]
    block_sizes = problem["block_sizes"]

    Y = numpy.zeros([num_constraints], float)
    line = f.readline()
    parts = line.split()
    for i in range(num_constraints):
        Y[i] = float(parts[i])

    Primal = []
    Dual = []

    # ignore Dual ############################ We work only with primal anyway

    for block_num in range(num_blocks):
        if block_sizes[block_num] < 0:
            Primal.append(matrix(RDF, 1, -block_sizes[block_num]))
            # Dual.append(numpy.zeros([1, -block_sizes[block_num]], float))
        else:
            Primal.append(matrix(RDF, block_sizes[block_num], block_sizes[block_num]))
            # Dual.append(numpy.zeros([block_sizes[block_num], block_sizes[block_num]], float))

    while True:
        line = f.readline()
        if len(line) == 0:
            break
        parts = line.split()
        if len(parts) == 0:
            continue
        constraint = int(parts[0])
        block_num = int(parts[1])
        row = int(parts[2])
        column = int(parts[3])
        value = float(parts[4])

        if constraint == 2:
            block = Primal[block_num - 1]
        else:
            continue
            # block = Dual[block_num - 1]

        if block_sizes[block_num - 1] < 0:
            block[0, column - 1] += value
        else:
            block[row - 1, column - 1] += value
            if row != column:
                block[column - 1, row - 1] += value

    solution = dict()
    solution["Y"] = Y
    solution["Primal"] = Primal
    # solution['Dual'] = Dual     # not loaded
    return solution


##################################################################### Rounding


def rational_approximation(V, den=roundingDenominator):
    """Entry-wise rounding of a matrix.
    Rounds every entry of V and multiplies by the den.
    Args:
        V: matrix to be rounded
        den: when rounding, every entry is multiplied by den
    Returns:
        matrix as 2D lists of rounded V, sclaled by den.
    """
    m = V.nrows()
    n = V.ncols()

    # W = matrix(field, m, n)
    W = [[field(0) for x in range(m)] for y in range(n)]
    for i in range(m):
        for j in range(n):
            # W[i,j] = QQ(float(real_part(V[i, j])))
            W[i][j] = field(round(real_part(V[i, j]) * den))
            # W[i,j] = field(round(real_part(V[i, j])*roundingDenominator))
    return W


# rounds one matrix A
def rounding_matrix(A):
    """Rounds a positive semidefinite matrix.
    Performs decomposition into eigenvalues and eigenvectors, they are
    simply rounded and mutliplied togethers. This assures A is positive
    semidefinite.
    Args:
        A: A square positive semidefinite matrix to be rounded
    Returns:
        Rounded matrix A as 2D list of integers with denominator
        roundingDenominator^3
    """

    n = A.nrows()
    evA, eigA = A.eigenmatrix_left()
    # print "Decomposition done"
    # evA is a diagonal matrix, eigenvalues on the diagonal
    # eigA is a matrix of eigenvectors

    # unlikely but possible that we get a complex numbers
    # as eigenvalues so we get rid of them
    for i in range(n):
        evA[i, i] = real_part(evA[i, i])
        if evA[i, i] < eigVZero:
            evA[i, i] = 0

    B = rational_approximation(eigA, roundingDenominator)
    evAR = rational_approximation(evA, roundingDenominator)
    # print "Rounding of decomposition done"
    C = [[0 for x in range(n)] for y in range(n)]
    # C = matrix(field, m, n)
    for i in range(n):
        for j in range(n):
            C[i][j] = sum(B[k][i] * evAR[k][k] * B[k][j] for k in range(n))
    # C = B.transpose() * evA * B
    # print "Product reconstructed"
    return C


# rounds one matrix A
def rounding_matrix_using_Cholesky(A):
    """Rounds a positive semidefinite matrix.
    Performs cholesky decomposition, rounds it and makes the matrix A back.
    Args:
    A: A square positive definite matrix to be rounded
    Returns:
    Rounded matrix A as 2D list of integers with denominator
    roundingDenominator^3
    """

    n = A.nrows()

    CholA = A.cholesky()

    CholARA = rational_approximation(CholA, roundingDenominator)

    # Test if diagonal entries are > 0
    for i in range(n):
        if CholARA[i][i] <= 0:
            print("Diagonal entry at ", i, " is ", CholARA[i][i])
            return

    C = [[0 for x in range(n)] for y in range(n)]
    for i in range(n):
        for j in range(n):
            C[i][j] = sum(
                roundingDenominator * CholARA[i][k] * CholARA[j][k] for k in range(n)
            )

    return C


@parallel
def roundAndPickleBlock(
    file, problem, solnMatrix, i, pickleBlock=True, useCholesky=False
):
    """Rounds one block of the solution matrix and saves it to a file.
    May run un parallel. Generates file file.result.rounded.Integer.i.pkl
    with rounded block.
    Args:
        file: Filename of the dat-s.
        problem: problem OUTLINE from the dat-s
        solnMatrix: List of all solution matrices.
            Actually only solnMatrix[i] is used.
        i: solnMatrix[i] is the one to process
    Returns:
        nothing returned
    """
    print("Working on rounding matrix", i, "/", problem["num_blocks"] - 1)
    if problem["block_sizes"][i] > 0:
        if useCholesky == False:
            rounded = rounding_matrix(solnMatrix[i])
        else:
            rounded = rounding_matrix_using_Cholesky(solnMatrix[i])
    elif problem["block_sizes"][i] != -problem["num_constraints"] - 1:
        m = solnMatrix[i].nrows()
        n = solnMatrix[i].ncols()
        rounded = [[0 for y in range(n)] for x in range(m)]
        # rounded = matrix(field,matrix.zero(m,n))
        for x in range(m):
            for y in range(n):
                rounded[x][y] = field(
                    real_part(
                        solnMatrix[i][x, y]
                        * roundingDenominator
                        * roundingDenominator
                        * roundingDenominator
                    )
                )
                if rounded[x][y] < eigVZero:
                    rounded[x][y] = 0
    else:
        m = solnMatrix[i].nrows()
        n = solnMatrix[i].ncols()
        rounded = [[0 for y in range(n)] for x in range(m)]
        # rounded = matrix(field,matrix.zero(m,n))
    if pickleBlock == True:
        output = open(file + ".result.rounded.Integer." + str(i) + ".pkl", "wb")
        pickle.dump(rounded, output)
        output.close()
        print("Matrix", i, "/", problem["num_blocks"] - 1, " is rounded and pickled")
    else:
        print("Matrix", i, "/", problem["num_blocks"] - 1, " is rounded")
    return rounded


def roundAndPickleBlocks_OLD(
    file, problem_outline, solnMatrix, useParallel=runParellel
):
    """Rounds and pickles all block in solnMatrix.
    Args:
        file: filename for dat-s file
        problem_outline: loaded header from file
        solnMatrix: blocks in solution
        useParallel: enable run in parallel
    Returns:
        nothing
    """
    if useParallel:
        list(
            roundAndPickleBlock(
                [
                    (file, problem_outline, solnMatrix, i)
                    for i in range(problem_outline["num_blocks"])
                ]
            )
        )
    else:
        for i in range(problem_outline["num_blocks"]):
            roundAndPickleBlock(file, problem_outline, solnMatrix, i)


def roundAndPickleBlocks(file, problem_outline, solnMatrix, useParallel=runParellel):
    """Rounds and pickles all block in solnMatrix.
    Args:
        file: filename for dat-s file
        problem_outline: loaded header from file
        solnMatrix: blocks in solution
        useParallel: enable run in parallel
    Returns:
        nothing
    """

    if useParallel:
        # In parallel processnig, the return value is a pair (input, output) instead of just output
        rounded_blocks_big = list(
            roundAndPickleBlock(
                [
                    (file, problem_outline, solnMatrix, i, False)
                    for i in range(problem_outline["num_blocks"])
                ]
            )
        )
        rounded_blocks = [x for tmp, x in rounded_blocks_big]
    else:
        rounded_blocks = list(
            [
                roundAndPickleBlock(file, problem_outline, solnMatrix, i, False)
                for i in range(problem_outline["num_blocks"])
            ]
        )

    output = open(file + ".result.rounded.Integer." + "all" + ".pkl", "wb")
    pickle.dump(rounded_blocks, output)
    output.close()
    print("Solution pickled.")


def unpickleRoundedMatrices(file, problem_outline):
    """Unpickles blocks pickled by function roundAndPickleBlocks.
    Args:
        file: filename for dat-s file
        problem_outline: Problem outline for file
    Returns:
        Array of unpickled blocks
    """
    input = open(file + ".result.rounded.Integer." + "all" + ".pkl", "rb")
    if not input.closed:
        roundedMatrices = pickle.load(input)
        input.close()
        return roundedMatrices

    roundedMatrices = []
    for i in range(problem_outline["num_blocks"]):
        input = open(file + ".result.rounded.Integer." + str(i) + ".pkl", "rb")
        roundedMatrices.append(pickle.load(input))
        input.close()
        print("Matrix", i + 1, "/", problem_outline["num_blocks"], "is unpickled")
    return roundedMatrices


##################################################### Testing rounded matrices


@parallel
def testOneConstraint(problem, solutionMatrices, scale, constants_scaled, i):
    """Tests one constraint in of the rounding.
    Can be executed in parallel.
    Args:
        problem: Loaded problem
        solutionMatrices: list of 2D arrays of solution
        scale: what is the scale for the objective of problem
        constants_scaled: Array of scaled constants from problem (all Integers)
        i: ID of constraint from problem to be tested
    Returns:
        pair [i,c] where c is evalutaion of the
    """

    num_constraints = problem["num_constraints"]
    num_blocks = problem["num_blocks"]
    block_sizes = problem["block_sizes"]
    blocks = problem["blocks"]
    # print "Have",num_constraints
    s = field(0)
    for j in range(0, num_blocks):
        if block_sizes[j] != 0 and block_sizes[j] != -num_constraints - 1:
            A = blocks[j][i + 1]
            x_len = len(solutionMatrices[j])
            y_len = len(solutionMatrices[j][0])
            # x_len = solutionMatrices[j].nrows()
            # y_len = solutionMatrices[j].ncols()
            if A.nrows() != x_len or A.ncols() != y_len:
                print("ERROR Solution:", x_len, y_len)
                print("A:", A.nrows(), A.ncols())
                continue
            # print "x_len=",x_len, "  y_len=",y_len
            # coordinates swapped for unit diagonal ones
            delta_s = sum(
                sum(A[x, y] * solutionMatrices[j][x][y] for y in range(y_len))
                for x in range(x_len)
            )
            s += delta_s
            # s += sum(sum(A.elementwise_product(solutionMatrices[j])))
            # print "Trying block",j,block_sizes[j],"  c-s=",(float(float(constants_scaled[i] -s)/(roundingDenominator^3))/scale),  "delta_s=",float((delta_s/(roundingDenominator^3))/scale)
    #        else:
    #            print "Trying block",j,block_sizes[j],"  c-s=",(float(float(constants_scaled[i] -s)/(roundingDenominator^3))/scale),  "delta_s=skipped"
    c = constants_scaled[i] - s
    print(
        i,
        "/",
        num_constraints,
        "rounded to",
        float(float(c) / (roundingDenominator ^ 3)) / scale,
    )
    return [i, c]


def testMatrices(
    problem,
    solutionMatrices,
    scale,
    test_from=0,
    test_to=10000000,
    useParellel=runParellel,
):
    """Tests solutionMatrices in constraints of a problem.
    Args:
        problem: loaded problem with constraints to test.
        solutionMatrices: solution in form of integer 2D array
        scale: Scale for the objective function
        test_from: First constraint to test
        test_to: First constrraint to not test
        useParallel: Allow run in parallel
    Returns:
        Scaled minimum found in the constraints.
    """
    loaded_constraints = problem["loaded_constraints"]
    num_constraints = problem["num_constraints"]
    num_blocks = problem["num_blocks"]
    block_sizes = problem["block_sizes"]
    blocks = problem["blocks"]
    constants = problem["constants"]
    min = 0
    min_i = 0
    min_real = 0
    const_scale = Integer(roundingDenominator) ^ 3
    constants_scaled = [Integer(c) * const_scale for c in constants]
    if test_to > loaded_constraints:
        test_to = loaded_constraints
    # print num_constraints
    if useParellel:
        for X, Y in list(
            testOneConstraint(
                [
                    (problem, solutionMatrices, scale, constants_scaled, i)
                    for i in range(test_from, test_to)
                ]
            )
        ):
            i = Y[0]
            c = Y[1]
            if min == 0 or c < min:
                min = c
                min_i = i
            min_real = float(float(min) / (roundingDenominator ^ 3))
    else:
        for i in range(test_from, test_to):
            x, c = testOneConstraint(
                problem, solutionMatrices, scale, constants_scaled, i
            )
            if min == 0 or c < min:
                min = c
                min_i = i
            min_real = float(float(min) / (roundingDenominator ^ 3))
    print(
        "Formal lower bound on (1 - obj. function) found at",
        min_i,
        "/",
        loaded_constraints,
        "and it is",
        float(min_real / scale),
    )
    return min


############################################################################ Testing solution from CSDP


def matrix_to_arrays(M):
    """Converst one matrix from solution to form expected by test.
    Args:
        M: RDF matrix
    Returns:
        2D array that is properly scaled as would come from roudning.
    """
    A = [
        [(roundingDenominator ^ 3) * M[x, y] for y in range(M.ncols())]
        for x in range(M.nrows())
    ]
    return A


def soln_matrix_to_arrays(solnMatrix):
    """Converts matrices loaded from solution to matrices that can be used in testing.
    Pretty much provides the required scale and form for matrices from solution
    to be tested the same way rounded matrices are.
    Args:
        solnMatrix: Array of solution matrices
    Returns:
        Array af 2D arrays created from solnMatrix
    """
    return [matrix_to_arrays(M) for M in solnMatrix]


def test_solution(file, scale, sdp_scaled=True, test_from=0, test_to=1000000):
    """Test matrices that were created by CSDP the same way rounded matrices are.
    Args:
        file: Filename with dat-s.
        scale: Scale of the objective for dat-s file.
        sdp_scaled: True if dat-s was already scaled before running CSDP.
        test_from: First constraint from program to test.
        test_to: First not tested. Testing is [test_from, test_to)
    Returns:
        Nothing. Prints the result of test.
    """
    print("Reading problem file...", file)
    if sdp_scaled:
        problem = read_problem_file(file, 1, test_from, test_to)
    else:
        problem = read_problem_file(file, scale, test_from, test_to)

    print("Reading solution file...")
    solution = read_solution_file(file + ".result", problem)

    print("Solution loaded....")
    solnMatrix = solution["Primal"]

    solnArray = soln_matrix_to_arrays(solnMatrix)

    print("Testing of rounded starts...")
    r = testMatrices(problem, solnArray, scale)

    print(
        r,
        "/",
        scale * (roundingDenominator ^ 3),
        " - ",
        float(r) / float(scale * (roundingDenominator ^ 3)),
    )


########################################################################## For easy programs


def only_rounding(file, scale, sdp_scaled=True):
    """Performs rounding and pickling of a solution from CSDP.
    Assume that the solution to dat-s file is called file.result
    Args:
        file: filename for dat-s file
        scale: What was the scale of the objective function to make
            it integral coefficients in dat-s file
        sdp_scaled: True if the dat-s file was already scaled.
    Returns:
        Nothing. Files with pickled rounded arrays are written to disk.
    """
    print("Reading problem outline...", file)
    print("With scale", scale)
    problem_outline = read_problem_outline(file)

    print("Reading solution file...")
    solution = read_solution_file(file + ".result", problem_outline)

    print("Solution loaded....")
    solnMatrix = solution["Primal"]

    print("Rounding using pickle starts...")
    roundAndPickleBlocks(file, problem_outline, solnMatrix)


def test_rounding(file, scale, sdp_scaled=True, test_from=0, test_to=10000000):
    """Unpickles and tests rounding"""
    print("Reading problem outline...", file)
    print("With scale", scale)
    problem_outline = read_problem_outline(file)

    print("Unpickling rounded matrices")
    roundedMatrices = unpickleRoundedMatrices(file, problem_outline)

    print("Reading problem file...", file)
    if sdp_scaled:
        problem = read_problem_file(file, 1, test_from, test_to)
    else:
        problem = read_problem_file(file, scale, test_from, test_to)

    print("Testing of rounded starts...")
    r = testMatrices(problem, roundedMatrices, scale)

    denom = scale * (roundingDenominator ^ 3)
    print(
        "Final upper bound on objective function is",
        -r + denom,
        "/",
        denom,
        "# = ",
        float(-r + denom) / float(denom),
    )
    return scale * (roundingDenominator ^ 3)


def round_program(file, scale, sdp_scaled=True):
    """Rounds and tests rounding.
    Args:
        file: Name of dat-s file to test
        scale: How much is the objective function scaled to get integer coefficients
        sdp_scaled: True if it was scaled before CSDP
    """
    # Ensure problem can actually be read before performing expensive computations.
    _ = read_problem_file(file, scale, 0, 1000000000)
    only_rounding(file, scale, sdp_scaled)
    test_rounding(file, scale, sdp_scaled)


def max_diff(arrayA, arrayB):
    maxdiff = 0
    for x in range(len(arrayA)):
        for y in range(len(arrayB)):
            if abs(arrayA[x][y] - arrayB[x][y]) > maxdiff:
                maxdiff = abs(arrayA[x][y] - arrayB[x][y])
    return maxdiff


def arrays_to_matrix(arrayA):
    """Converst one array from rounding to a matrix that looks like from CSDP
    Args:
    arrayA: 2D array that is properly scaled as would come from roudning.
    Returns:
    RDF matrix
    """
    M = matrix(RDF, len(arrayA), len(arrayA[0]))
    for x in range(M.nrows()):
        for y in range(M.ncols()):
            M[x, y] = arrayA[x][y] / (roundingDenominator ^ 3)
    #    A = [[(roundingDenominator^3)*M[x,y] for y in range(M.ncols())] for x in range(M.nrows())]
    return M

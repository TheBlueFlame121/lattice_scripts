from tqdm import tqdm

########################
##### DEFINITIONS ######
########################
q = 2**23 - 2**13 + 1
Zq = Zmod(q)
Zqx = PolynomialRing(Zq, 'x')
x = Zqx.gen()
prim_gen = Zq(10)
n = 256
# Rq = Zqx.quotient(x^256 + 1)
# root_of_unity = Zq(1753)

def br(i, n):
    return int(bin(i)[2:].zfill(n)[::-1], 2)




##########################
##### Base Functions #####
##########################
def fft(P, wn):
    n = len(P)
    assert Integer(n).is_power_of(2)
    if n == 1:
        return P
    Pe = P[::2]
    ye = fft(Pe, wn**2)
    Po = P[1::2]
    yo = fft(Po, wn**2)

    y = [0]*n
    w = 1
    for j in range(n/2):
        y[j] = ye[j] + w*yo[j]
        y[j + n/2] = ye[j] - w*yo[j]
        w*=wn
    return y


def ifft(P, wn):
    y = fft(P, 1/wn)
    return [i/len(P) for i in y]


def fft_iter(P):
    n = len(P)
    assert Integer(n).is_power_of(2)

    w = pow(prim_gen, (q-1)/(2*n))
    zetas = [pow(w, br(i, log(n, 2))) for i in range(n)]

    PairsInGroup = n/2
    NumOfGroups = 1
    Distance = n/2
    while NumOfGroups < n:
        for k in range(NumOfGroups):
            JFirst = 2*k*PairsInGroup
            JLast = JFirst + PairsInGroup - 1
            Jtwiddle = k
            W = zetas[Jtwiddle]
            for j in range(JFirst, JLast+1):
                temp = W*P[j+Distance]
                P[j+Distance] = P[j] - temp
                P[j] = P[j] + temp
        PairsInGroup/=2
        NumOfGroups*=2
        Distance/=2
    return P


def ifft_iter(P):
    n = len(P)
    assert Integer(n).is_power_of(2)

    w = pow(prim_gen, (1-q)/(2*n))
    inv_zetas = [pow(w, br(i, log(n, 2))) for i in range(n)]

    NumOfProblems = 1
    ProblemSize = n
    Distance = 1
    while ProblemSize > 1:
        for JFirst in range(NumOfProblems):
            J = JFirst
            Jtwiddle = 0
            while J < n-1:
                W = inv_zetas[Jtwiddle]
                Temp = P[J]
                P[J] = Temp + P[J+Distance]
                P[J+Distance] = (Temp-P[J+Distance])*W
                Jtwiddle+=1
                J += 2*NumOfProblems
        NumOfProblems*=2
        ProblemSize/=2
        Distance*=2
    return [i/n for i in P]


def ntt(P):
    n = len(P)
    psi = pow(prim_gen, (q-1)/(2*n))

    tP = [(psi**i)*j for i,j in enumerate(P)]
    return fft(tP, psi**2)


def intt(P):
    n = len(P)
    psi_inv = pow(prim_gen, (1-q)/(2*n))

    tP = ifft(P, 1/psi_inv**2)
    return [(psi_inv**i)*j for i, j in enumerate(tP)]


def ntt_iter(P):
    n = len(P)
    psi = pow(prim_gen, (q-1)/(2*n))

    tP = [(psi**i)*j for i,j in enumerate(P)]
    return fft_iter(tP)


def intt_iter(P):
    n = len(P)
    psi_inv = pow(prim_gen, (1-q)/(2*n))

    tP = ifft_iter(P)
    return [(psi_inv**i)*j for i, j in enumerate(tP)]




##########################
##### Test Functions #####
##########################
def test_fft(degree, number):
    root_of_unity = pow(prim_gen, (q-1)/degree)
    for _ in tqdm(range(number)):
        P = Zqx.random_element(degree-1)
        assert [P(root_of_unity**i) for i in range(degree)] == fft(P.list(), root_of_unity)

def test_ifft(degree, number):
    root_of_unity = pow(prim_gen, (q-1)/degree)
    for _ in tqdm(range(number)):
        P = Zqx.random_element(degree-1)
        assert ifft(fft(P.list(), root_of_unity), root_of_unity) == P.list()

def test_fft_iter(degree, number):
    root_of_unity = pow(prim_gen, (q-1)/degree)
    for _ in tqdm(range(number)):
        P = Zqx.random_element(degree-1)
        assert set(fft(P.list(), root_of_unity)) == set(fft_iter(P.list()))

def test_ifft_iter(degree, number):
    for _ in tqdm(range(number)):
        P = Zqx.random_element(degree-1)
        assert ifft_iter(fft_iter(P.list().copy())) == P.list()

def test_ntt(degree, number):
    root_of_unity = pow(prim_gen, (q-1)/(2*degree))
    for _ in tqdm(range(number)):
        P = Zqx.random_element(degree-1)
        assert [P(root_of_unity**(2*i + 1)) for i in range(degree)] == ntt(P.list())

def test_intt(degree, number):
    for _ in tqdm(range(number)):
        P = Zqx.random_element(degree-1)
        assert intt(ntt(P.list())) == P.list()

def test_ntt_iter(degree, number):
    root_of_unity = pow(prim_gen, (q-1)/(2*degree))
    for _ in tqdm(range(number)):
        P = Zqx.random_element(degree-1)
        assert set(P(root_of_unity**(2*i + 1)) for i in range(degree)) == set(ntt_iter(P.list()))

def test_intt_iter(degree, number):
    for _ in tqdm(range(number)):
        P = Zqx.random_element(degree-1)
        assert intt_iter(ntt_iter(P.list())) == P.list()




####################################
##### Multiplication Functions #####
####################################
def normal_convolution(degree):
    P = Zqx.random_element(degree-1)
    Q = Zqx.random_element(degree-1)
    R = (P*Q) % (x^degree - 1)


def fft_convolution(degree):
    root_of_unity = pow(prim_gen, (q-1)/degree)
    P = fft(Zqx.random_element(degree-1).list(), root_of_unity)
    Q = fft(Zqx.random_element(degree-1).list(), root_of_unity)
    R = ifft([i*j for i, j in zip(P, Q)], root_of_unity)


def fft_iter_convolution(degree):
    P = fft_iter(Zqx.random_element(degree-1).list())
    Q = fft_iter(Zqx.random_element(degree-1).list())
    R = ifft_iter([i*j for i, j in zip(P, Q)])


def normal_mul(degree):
    P = Zqx.random_element(degree-1)
    Q = Zqx.random_element(degree-1)
    R = (P*Q) % (x^degree + 1)


def ntt_mul(degree):
    P = ntt(Zqx.random_element(degree-1).list())
    Q = ntt(Zqx.random_element(degree-1).list())
    R = Zqx(intt([i*j for i, j in zip(P, Q)]))

def ntt_iter_mul(degree):
    P = ntt_iter(Zqx.random_element(degree-1).list())
    Q = ntt_iter(Zqx.random_element(degree-1).list())
    R = Zqx(intt_iter([i*j for i, j in zip(P, Q)]))

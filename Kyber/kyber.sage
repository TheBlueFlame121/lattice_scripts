from tqdm import tqdm

# q = 7681
# n = 256
q = 17
n = 4

Zqx.<x> = Zmod(q)[]
R.<x> = Zqx.quotient(x^n + 1)

def set_vals(_q, _n):
    global q, n, Zqx, R
    q = _q
    n = _n
    Zqx.<x> = Zmod(q)[]
    R.<x> = Zqx.quotient(x^n + 1)

def ternary_polynomial():
    a = R.random_element()
    return R([ZZ(i)%3 - 1 for i in a.list()])

def binary_polynomial():
    a = R.random_element()
    return R([ZZ(i)%2 for i in a.list()])


def generate_keys(k:int = 2):
    A_list = [R.random_element() for i in range(k*k)]
    A = Matrix(k, k, A_list)
    
    s_list = [ternary_polynomial() for i in range(k)]
    s = vector(s_list)

    e_list = [ternary_polynomial() for i in range(k)]
    e = vector(e_list)

    t = A*s + e
    return s, (A, t)


def encrypt(m, pub):
    A, t = pub
    k = t.length()

    r = vector([ternary_polynomial() for i in range(k)])
    e1 = vector([ternary_polynomial() for i in range(k)])
    e2 = ternary_polynomial()

    m*=round(q/2)

    u = A.transpose()*r + e1
    v = t*r + e2 + m

    return (u, v)


def decrypt(cip, priv):
    u, v = cip
    s = priv

    m_n = v - s*u
    m = []
    for i in map(int, m_n.list()):
        if abs(i - round(q/2))<=min(abs(i), abs(i-q)):
            m.append(round(q/2))
        else:
            m.append(0)

    # m = [round(q/2) if abs(i - round(q/2))<=min(abs(i), abs(i-q)) else 0 for i in map(int, m_n.list())]

    return R(m)/round(q/2)


def test_correctness():
    priv, pub = generate_keys()
    m = binary_polynomial()
    cip = encrypt(m, pub)
    m_prime = decrypt(cip, priv)
    return m == m_prime

def test_percentage(times):
    count = 0
    for i in tqdm(range(times)):
        if test_correctness():
            count += 1
    return count

# TODO:
# Code overview
# Correctness
# Efficiency improvements using symmetric constructs
# Implementation benefits/Parameter sets


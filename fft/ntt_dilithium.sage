from tqdm import tqdm

q = 2**23 - 2**13 + 1
Zq = Zmod(q)
Zqx = PolynomialRing(Zq, 'x')
x = Zqx.gen()
n = 256
n_inv = Zq(1/n)
psi = Zq(1753)

def br(i, n):
    return int(bin(i)[2:].zfill(n)[::-1], 2)

zetas = [pow(psi, br(i, log(n, 2))) for i in range(n)]
inv_zetas = [pow(psi, -br(i, log(n, 2))) for i in range(n)]

def Cooley_Tukey_NTT(x):
    t = 128
    m = 1
    while m < 256:
        k = 0
        for i in range(m):
            w = zetas[m + i]
            for j in range(k, k+t):
                v = x[j + t]*w
                x[j + t] = x[j] - v
                x[j] = x[j] + v
            k += 2*t
        t /= 2
        m *= 2
    return x

def Gentlemen_Sande_INTT(x):
    t = 1
    m = 128
    while m > 0:
        k = 0
        for i in range(m):
            S = inv_zetas[m + i]
            for j in range(k, k+t):
                U = x[j]
                V = x[j + t]
                W = U - V
                x[j] = U + V
                x[j + t] = W*S
            k += 2*t
        t *= 2
        m //= 2
    for i in range(256):
        x[i] *= n_inv
    return x


def ntt(coeffs):
    k, l = 0, 128
    while l > 0:
        start = 0
        while start < 256:
            k = k + 1
            zeta = zetas[k]
            for j in range(start, start + l):
                t = zeta*coeffs[j+l]
                coeffs[j+l] = coeffs[j] - t
                coeffs[j]   = coeffs[j] + t
            start = l + (j + 1)
        l >>= 1
    return coeffs

def intt(coeffs):
    l, k = 1, 256
    while l < 256:
        start = 0
        while start < 256:
            k = k - 1
            zeta = -zetas[k]
            for j in range(start, start+l):
                t = coeffs[j]
                coeffs[j]   = t + coeffs[j+l]
                coeffs[j+l] = t - coeffs[j+l]
                coeffs[j+l] = zeta*coeffs[j+l]
            start = j + l + 1
        l = l << 1
    for j in range(256):
        coeffs[j] = coeffs[j]/n
    return coeffs


def verify_textbook():
    P = Zqx.random_element(255)
    Q = Zqx.random_element(255)
    R = (P*Q)%(x^256+1)
    p = Cooley_Tukey_NTT(P.list())
    q = Cooley_Tukey_NTT(Q.list())
    r = Gentlemen_Sande_INTT([i*j for i, j in zip(p, q)])
    assert r == R.list()

def verify_ref():
    P = Zqx.random_element(255)
    Q = Zqx.random_element(255)
    R = (P*Q)%(x^256+1)
    p = ntt(P.list())
    q = ntt(Q.list())
    r = intt([i*j for i, j in zip(p, q)])
    assert r == R.list()

def test():
    p = ntt(Zqx.random_element(255).list())
    q = ntt(Zqx.random_element(255).list())
    r = Zqx(intt([i*j for i, j in zip(p, q)]))

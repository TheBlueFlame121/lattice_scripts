from tqdm import tqdm

q = 3329
Zq = Zmod(q)
Zqx = PolynomialRing(Zq, 'x')
x = Zqx.gen()
n = 256
n_inv = Zq(1/n)
psi = Zq(17)

def br(i, n):
    return int(bin(i)[2:].zfill(n)[::-1], 2)

zetas = [pow(psi, br(i, 7)) for i in range(128)]

def ntt(coeffs):
    k, l = 1, 128
    while l >= 2:
        start = 0
        while start < 256:
            zeta = zetas[k]
            k = k + 1
            for j in range(start, start + l):
                t = zeta*coeffs[j+l]
                coeffs[j+l] = coeffs[j] - t
                coeffs[j]   = coeffs[j] + t
            start = l + (j + 1)
        l >>= 1
    return coeffs

def intt(coeffs):
    l, k = 2, 127
    while l <= 128:
        start = 0
        while start < 256:
            zeta = -zetas[k]
            k = k - 1
            for j in range(start, start+l):
                t = coeffs[j]
                coeffs[j]   = t + coeffs[j+l]
                coeffs[j+l] = t - coeffs[j+l]
                coeffs[j+l] = zeta*coeffs[j+l]
            start = j + l + 1
        l = l << 1
    for j in range(256):
        coeffs[j] = coeffs[j]/128
    return coeffs

def basemul(a, b, zeta):
    r0 = (a[1]*b[1])*zeta + a[0]*b[0]
    r1 = a[0]*b[1] + b[0]*a[1]
    return [r0, r1]

def poly_basemul(a, b):
    res = []
    for i in range(64):
        r0, r1 = basemul(
                a[4*i: 4*i+2],
                b[4*i: 4*i+2],
                zetas[64+i])
        r2, r3 = basemul(
                a[4*i+2: 4*i+4],
                b[4*i+2:4*i+4],
                -zetas[64+i])
        res += [r0, r1, r2, r3]
    return res

def verify_ref():
    P = Zqx.random_element(255)
    Q = Zqx.random_element(255)
    R = (P*Q)%(x^256+1)
    p = ntt(P.list())
    q = ntt(Q.list())
    r = intt(poly_basemul(p, q))
    assert r == R.list()

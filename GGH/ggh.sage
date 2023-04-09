import random

def hadamard_ratio(basis):
    n = basis.dimensions()[0]
    if not basis.is_square():
        raise Exception("no non-square matrices allowed")
    d = abs(basis.determinant())
    l = prod([v.norm() for v in basis])
    return RR((d/l)**(1/n))


def random_unimodular_matrix(n, d):
    up = random_matrix(ZZ, n, x=-d, y=d+1).list()
    low = random_matrix(ZZ, n, x=-d, y=d+1).list()
    for i in range(n):
        for j in range(n):
            if i < j :
                up[i*n + j] = 0
                continue
            if i == j :
                up[i*n + j] = random.choice([-1, 1])
                low[i*n + j] = random.choice([-1, 1])
                continue
            if i > j :
                low[i*n + j] = 0
                continue
    return Matrix(ZZ, n, n, up)*Matrix(ZZ, n, n, low)
    


def generate_private(n, d):
    while True:
        priv = d * identity_matrix(n)
        priv[0] = vector(ZZ, [randint(-d, d) for _ in range(n)])
        if priv.is_singular() or hadamard_ratio(priv) < 0.9:
            continue
        return priv


def generate_public(priv):
    n = priv.dimensions()[0]
    public = random_unimodular_matrix(n, 10) * priv
    while hadamard_ratio(public) > 0.01:
        public = random_unimodular_matrix(n, 10) * public
    return public


def get_message(n, d):
    return random_vector(ZZ, n, x=-d, y=d+1)


def encrypt(pub, m, e):
    n = pub.dimensions()[0]
    r = random_vector(ZZ, n, x=-e, y=e+1)
    return (m*pub) + r


def decrypt(c, priv, pub):
    mw = cvp_rounding(c, priv)
    return mw*pub.inverse()


def cvp_rounding(target, basis):
    t = target.change_ring(RR)
    br = basis.change_ring(RR)
    c = vector(ZZ, [round(i) for i in t*br.inverse()])
    return c*basis


def test_correctness(n, d, e):
    a = generate_private(n, d)
    b = generate_public(a)
    m = get_message(n, d)
    c = encrypt(b, m, e)
    return decrypt(c, a, b) == m


def test_security(n, d, e):
    a = generate_private(n, d)
    b = generate_public(a)
    m = get_message(n, d)
    c = encrypt(b, m, e)
    return decrypt(c, b, b) != m

# To do:
# 1. Code overview
# 2. Parameter generation
# 3. correctess and security test
# 4. Cvp effectiveness by comparing distance
# 5. Variation in error

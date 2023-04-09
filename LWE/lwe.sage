def generate_keys(n, m, q):
    A_list = [randint(0, q-1) for i in range(n*m)]
    A = Matrix(Zmod(q), m, n, A_list)

    s_list = [randint(0, q-1) for i in range(n)]
    s = vector(Zmod(q), s_list)

    e_list = [randint(-2, 2) for i in range(m)]
    e = vector(Zmod(q), e_list)

    b = A*s + e
    return s, (A, b)

def encrypt(message, pub, n, m, q):
    assert message in [0, 1] and "Message should be only one bit"

    A, b = pub
    x = vector([randint(-1, 1) for i in range(m)])

    message*=round(q/2)

    u = x*A
    v = x*b + message
    
    return (u, v)

def decrypt(cip, priv, n, q):
    u, v = cip
    s = priv
    m_n = int(v - s*u)

    if abs(m_n - round(q/2))<=min(abs(m_n), abs(m_n-q)):
        m = 1
    else:
        m = 0

    return m

def test_correctness():
    # Parameters
    n = 64
    q = next_prime(0x10001)
    m = 96

    priv, pub = generate_keys(n, m, q)
    mess = randint(0, 1)
    cip = encrypt(mess, pub, n, m, q)
    mess_prime = decrypt(cip, priv, n, q)
    return mess == mess_prime


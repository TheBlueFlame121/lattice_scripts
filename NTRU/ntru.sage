Zx.<x> = ZZ[]

def convolution(f, g, n):
      return (f * g) % (x^n-1)

def balanced_mod(f,q):
  g = list(((i + q//2) % q) - q//2 for i in f.list())
  return Zx(g)

def random_poly(n, d):
    assert d <= n
    result = n*[0]
    for j in range(d):
        while True:
            r = randrange(n)
            if not result[r]: break
        result[r] = 1-2*randrange(2)
    return Zx(result)

def poly_inv_prime(f, p, n):
    T = Zx.change_ring(Zmod(p)).quotient(x^n-1)
    return Zx(lift(1 / T(f)))

def poly_inv_pow2(f, q, n):
    assert q.is_power_of(2)
    g = poly_inv_prime(f, 2, n)
    while True:
        r = balanced_mod(convolution(g, f, n),q)
        if r == 1: return g
        g = balanced_mod(convolution(g, 2-r, n),q)

def keypair(n, p, q, d):
    while True:
        try:
            f = random_poly(n, d)
            fp = poly_inv_prime(f,p,n)
            fq = poly_inv_pow2(f,q,n)
            break
        except InterruptedError:
            exit()
        except:
            pass
    g = random_poly(n, d)
    publickey = balanced_mod(p * convolution(fq, g, n),q)
    secretkey = f,fp
    return publickey,secretkey

def random_message(n):
    result = list(randrange(3) - 1 for j in range(n))
    return Zx(result)

def encrypt(m, pub, n, q, d):
    r = random_poly(n, d)
    return balanced_mod(convolution(pub, r, n) + m, q)

def decrypt(c, priv, n, p, q):
    f,fp = priv
    a = balanced_mod(convolution(c, f, n), q)
    return balanced_mod(convolution(a,fp, n), p)


def test_correctness(n, p, q, d):
    pub, priv = keypair(n, p, q, d)
    m = random_message(n)
    c = encrypt(m, pub, n, q, d)
    print(m == decrypt(c, priv, n, 3, q))

for i in range(10):
    test_correctness(743, 3, 2048, 495)

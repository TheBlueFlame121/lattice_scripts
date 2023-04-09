def euclids(a, b):
    if b == 0:
        return a
    print(f"({a}) = ({a//b}) * ({b}) + ({a%b})")
    return euclids(b, a%b)

def ext_euclids(a, b):
    if b == 0:
        return 1, 0, a
    print(f"({a}) = ({a//b}) * ({b}) + ({a%b})")
    x1, y1, g = ext_euclids(b, a%b)
    x = y1
    y = x1 - y1*(a//b)
    print(f"({a})*({x}) + ({b})*({y}) = ({g}) ")
    return x, y, g

def mod_inv(a, b):
    x, y, g = ext_euclids(a, b)
    assert g == 1 and "Modular Inverse does not exist"
    return x

def euclids_poly(a, b):
    assert a.is_monic() and b.is_monic()
    print(f"({a}) = ({a//b}) * ({b}) + ({a%b})")
    if (a%b) == 0:
        return b
    return euclids_poly(b, (a%b).monic())

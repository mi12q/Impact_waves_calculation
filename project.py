import math
import numpy as np

def func(coef):
    return np.poly1d(coef)

def newtons_method(a,b,coef,eps):
    f = func(coef)
    df = f.deriv()
    x = (a+b)/2
    h = f(x)/df(x)
    while abs(h) > eps:
        h = f(x) / df(x)
        x = x - h
    return x

def furie_sequence(coef):
    f = func(coef)
    sequence = [f]
    n = len(coef)
    i = 0
    while i < n-1:
        f = f.deriv()
        sequence.append(f)
        i += 1
    return sequence

def sturm_sequence(coef):
    f = func(coef)
    sequence = [f, f.deriv()]
    n = len(coef)
    i = 2
    while i < n:
        p = np.polydiv(sequence[i-2],sequence[i-1])
        new_coeffs = p[1].coeffs * (-1)
        new_f = func(new_coeffs)
        sequence.append(new_f)
        i += 1
    return sequence

def sign(x):
    if x >= 0:
        return 1
    else:
        return -1

def variation(seq):
    n = 0
    for i in range(1,len(seq)):
        if seq[i-1] == seq[i] * (-1):
            n+=1
    return n

def var(sequence,x):
    V1 = []
    for polynom in sequence:
        V1.append(sign(polynom(x)))
    var1 = variation(V1)
    return var1

def calculate_interval(coef):
    abs_coef = []
    for i in coef:
        abs_coef.append(abs(i))
    A = max(abs_coef[1:])
    B = max(abs_coef[:-1])
    x1 = abs_coef[-1]/(abs_coef[-1]+B)
    x2 = 1 + A/abs(abs_coef[0])
    return [x1,x2]

def solve_equation(a,b,coef,eps):
    b0 = b
    sequence = sturm_sequence(coef)
    roots = []

    if abs(var(sequence,a)-var(sequence,b)) == 0:
        return None
    elif abs(var(sequence,a)-var(sequence,b)) == 1:
        return newtons_method(a, b, coef,10**(-6))

    elif abs(var(sequence,a)-var(sequence,b)) > 1:
        while abs(var(sequence,a)-var(sequence,b))!=1 :
            c = (a+b)/2
            while var(sequence, a) - var(sequence, c) == 0:
                    a = c
                    c = (a+b)/2
            while var(sequence, b) - var(sequence, c) == 0:
                    b = c
                    c = (a+b)/2
            while abs(var(sequence,a)-var(sequence,b)) != 1:
                b = (a+b)/2
            root = newtons_method(a,b,coef,eps)
            roots.append(root)
            a = b
            b = b0
    root = newtons_method(a, b, coef, eps)
    roots.append(root)

    return roots


# p = solve_equation(-4,4,[1,2,-3],10**(-16))
# print(p)
# p = solve_equation(-5,4,[1,5,5,-5,-6],10**(-16))
# print(p)
#print(sturm_sequence([1,1,-1,-1]))


import math as math
import project as pr

def calculate_vel(y0,y3,ro0,ro3,U0,U3,P0,P3):
    C0 = math.sqrt(y0*P0/ro0)
    C3 = math.sqrt(y3*P3/ro3)
    X = P3/P0
    # Y = P1/P0
    alpha0 = (y0 +1)/(y0-1)
    alpha3 = (y3+1)/(y3-1)
    e3 = 2*C3**2/(y3*(y3-1)*(U3-U0)**2)
    e0 = 2 * C0 ** 2 / (y0 * (y0 - 1) * (U3 - U0)**2)
    a0 = (alpha0*e3 - alpha3*X*e0)**2
    a1 = 2*((alpha0*e3-alpha3*X*e0)*(e3*(1-2*alpha0*X)-e0*X*(X-2*alpha3))-alpha0*alpha3*X*(alpha0*e3+alpha3*X*e0))
    a2 = (e3**2*(6*alpha0**2*X**2-8*alpha0*X+1)-2*e0*e3*X*(alpha0*alpha3*(X**2+4*X+1)-2*(X+1)*(alpha3+alpha0*X)+X)+e0**2*X**2*(6*alpha3**2-8*alpha3*X+X**2)
          +alpha0**2*alpha3**2*X**2-2*alpha0*X*e3*(alpha0*X-2*alpha0*alpha3*X+2*alpha3)-2*alpha3*X**2*e0*(alpha3+2*alpha0*X-2*alpha0*alpha3))

    a3 = -2*X*(2*e3**2*(alpha0**2*X**2-3*alpha0*X+1)+e0*e3*((alpha3+alpha0*X)*(X**2+4*X+1)-2*alpha0*alpha3*X*(X+1)-2*X*(X+1))+
              2*e0**2*X*(X**2-3*alpha3*X+alpha3**2)-alpha0*alpha3*X*(alpha0*X+alpha3)+e3*(alpha0**2*alpha3*X**2-2*X*(2*alpha0*alpha3+alpha0**2*X)+(2*alpha0*X+alpha3))
               +e0*X*(alpha0*alpha3**2-2*alpha3*(alpha3+2*alpha0*X)+2*alpha3*X+alpha0*X**2))
    a4 = X**2*(e3**2*(alpha0**2*X**2-8*alpha0*X+6)-2*e0*e3*(alpha0*alpha3*X-2*(X+1)*(alpha3+alpha0*X)+X**2+4*X+1)+
               e0**2*(alpha3**2-8*alpha3*X+6*X**2)+(alpha3**2+4*alpha0*alpha3*X+alpha0**2*X**2)-2*e3*((alpha0**2*X+2*alpha0*alpha3)*X-2*(2*alpha0*X+alpha3)+1)-
          2*e0*(alpha3*(2*alpha0*X+alpha3)-2*X*(2*alpha3+alpha0*X)+X**2))
    a5 = 2*X**3*(e3**2*(alpha0*X-2)-e0*e3*(alpha0*X-2+alpha3-2*X)+e0**2*(alpha3-2*X)+(alpha3+alpha0*X)-e3*(2*alpha0*X+alpha3-2)-e0*(2*alpha3+alpha0*X-2*X))
    a6 = X**4*((e3-e0)**2+1-2*(e3+e0))
    coeffs = [a0,a1,a2,a3,a4,a5,a6]
    print('Коефициенты:', ', '.join(map(str, coeffs)))
    print("\n")
    interval = pr.calculate_interval(coeffs)
    Y = pr.solve_equation(interval[0],interval[1],coeffs,10**(-10))
    D = []

    for y in Y:
        p1 = y*P0
        u1 = U0 + (C0 * math.sqrt(2)/(math.sqrt(y0 * (y0 - 1)))) * (y - 1) / math.sqrt(1 + alpha0 * y)
        d3 = U3 - ((p1-P3)/ro3)*1/(U3-u1)
        d0 = U0 - ((p1 - P0) / ro0) * 1 / (U0 - u1)
        D.append((d0, d3))
        print("Root:",y )
        print("U1:", u1)
        print("D0:", d0)
        print("D3:", d3)

        u1 = U0 - (C0 * math.sqrt(2)/(math.sqrt(y0 * (y0 - 1)))) * (y - 1) / math.sqrt(1 + alpha0 * y)
        d3 = U3 - ((p1-P3)/ro3)*1/(U3-u1)
        d0 = U0 - ((p1 - P0) / ro0) * 1 / (U0 - u1)
        D.append((d0, d3))
        print("U1:", u1)
        print("D0:", d0)
        print("D3:", d3)
        print("\n")

    return D

calculate_vel(5/3,5/3,11.37,7.9,-2.28*10**(4),2.72*(10)**4,1.17928*10**(9),3.04*10**(9))
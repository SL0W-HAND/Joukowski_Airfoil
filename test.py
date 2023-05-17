import numpy as np
gamma = 1
eta0 = complex(0.2, 0.15)
u0 = 100
gamma = 4*np.pi*u0*(eta0.imag)
a = 1.2
R = 1


def potential(x, y):
    res = gamma*np.log(np.abs(eta0 - x/2 - complex(0, 1)*y/2 + np.sqrt(-4*a**2 + x**2 + 2*complex(0, 1)*x*y - y**2)/2))/(2*np.pi) + R**2*u0*(-y/2 + (4*x**2*y**2 + (-4*a**2 + x**2 - y**2)**2)**(1/4)*np.sin(np.arctan2(2*x*y, -4*a**2 + x**2 - y**2)/2)/2 + (eta0.imag))/((x/2 - (4*x**2*y**2 + (-4*a**2 + x**2 - y**2)**2)**(1/4)*np.cos(
        np.arctan2(2*x*y, -4*a**2 + x**2 - y**2)/2)/2 - (eta0.real))**2 + (y/2 - (4*x**2*y**2 + (-4*a**2 + x**2 - y**2)**2)**(1/4)*np.sin(np.arctan2(2*x*y, -4*a**2 + x**2 - y**2)/2)/2 - (eta0.imag))**2) + u0*(y/2 - (4*x**2*y**2 + (-4*a**2 + x**2 - y**2)**2)**(1/4)*np.sin(np.arctan2(2*x*y, -4*a**2 + x**2 - y**2)/2)/2 - (eta0.imag))
    return res


def derivate(func, x0, y0, param):
    h = 0.01
    if param == "x":
        print("x")
        return (func(x0+h, y0) - func(x0, y0))/h
    else:
        return (func(x0, y0+h) - func(x0, y0))/h


print(derivate(potential, 1, 1, "x"))

import numpy as np
from sympy import *
x, y, u0, r, Gamma, a, pi, zeta_x, zeta_y = symbols(
    'x y u0 R Gamma a  pi zeta_x  zeta_y', real=True)
z, eta, zeta, eta0 = symbols('z eta zeta eta0')

#z = x+I*y

zeta = x+I*y

z = (zeta-sqrt(zeta**2-4*a**2))/2-eta0


potential = u0*z+(u0*r**2)/z + (I*Gamma*log(z))/(2*pi)


# print(im(potential))

#print('stream function')
# print(latex(im(potential)))

# print('vx')
#print(diff(re(potential), x))
# print('vy')
#print(diff(re(potential), y))

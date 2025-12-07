import cmath
import numpy as np
# Define the parameters of the Joukowski airfoil
a = 1.2  # radius of the circle
mu = complex(0, 0)  # center of the circle

# Define the complex point to be checked
p = complex(0, 0)
p = -mu
# Calculate the Joukowski transform of the airfoil


def inverseShiftfunction(eta):
    return eta - mu


def joukowski(z):
    return z+(a**2)/z

# Calculate the transformed coordinates of the complex point


# Check if the magnitude of the transformed complex point is less than or equal to 1
if np.sqrt((((inverseShiftfunction(p))).real)**2 + (inverseShiftfunction(p).imag)**2) < 1:
    print("The point is inside the Joukowski airfoil.")
else:
    print("The point is outside the Joukowski airfoil.")

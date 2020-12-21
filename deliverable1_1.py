import sys
import numpy as np
from ase import Atoms
from ase.units import eV, Ang, GPa
import matplotlib.pyplot as plt

sys.path.append(r'files')
import Morse

# calculator object
calc = Morse.MorsePotential()

# range chosen to show negative region of morse_potential against distance
distance = np.linspace(2, 5, 50)

#initialise lists
potential_energy = []
forces = []

# calculate morse potential
for dist in distance:
    d = dist * Ang
    a = Atoms('2Cu', positions=[(0., 0., 0.), (0., 0., d)])

    a.set_calculator(calc)
    potential_energy.append(a.get_potential_energy())

plt.plot(distance, potential_energy)
plt.xlabel('Distance/ Angstroms')
plt.ylabel('Morse potential/ eV')
plt.title("Morse potential against distance")
plt.show()

# calculates magnitude of force between two atoms
for dist in distance:
    d = dist * Ang
    a = Atoms('2Cu', positions=[(0., 0., 0.), (0., 0., d)])

    a.set_calculator(calc)
    force = a.get_forces()

    # magnitude of the force on one atom
    forces.append(np.linalg.norm(force[0]))

plt.plot(distance, forces)
plt.xlabel('Distance/ Angstroms')
plt.ylabel('Force')
plt.title("Force against distance")
plt.show()

# unit tests...

# forces as obtained from get_forces()
forces_test = []
gradient_test = []

# get_force method
for dist in distance:
    d = dist * Ang
    a = Atoms('2Cu', positions=[(0., 0., 0.), (0., 0., d)])

    a.set_calculator(calc)
    force = a.get_forces()

    # force on one atom (chose this one so graphs can be compared easily)
    forces_test.append(force[0][2])

# gradient of Pot. energy method
for dist in distance:
    epsilon = 1e-9
    # pot2 is pot energy at position (current + epsilon d); pot1 is at current position
    d = dist * Ang
    a2 = Atoms('2Cu', positions=[(0., 0., 0.), (0., 0., d + epsilon)])
    a1 = Atoms('2Cu', positions=[(0., 0., 0.), (0., 0., d)])
    a2.set_calculator(calc)
    a1.set_calculator(calc)
    pot2 = a2.get_potential_energy()
    pot1 = a1.get_potential_energy()

    # gradient of the potential energy function
    tempGrad = (pot2 - pot1) / epsilon
    gradient_test.append(tempGrad)

# comparison graphs
plt.plot(distance, forces_test)
plt.plot(distance, gradient_test, '--')
plt.xlabel('Distance/ Angstroms')
plt.ylabel('Force/ eV')
plt.title("Force vs distance")
plt.show()

# checks actual values - prints the largest difference between the two
forceDiffArray = abs(np.array(forces_test) - np.array(gradient_test))
for i in forceDiffArray:
    assert (i < 1 * 10 ** -3)
import sys
import numpy as np
import matplotlib.pyplot as plt
from ase.build import bulk
from ase.units import GPa

sys.path.append(r'files')
import Morse

calc = Morse.MorsePotential()

epsilon = 0.02
strains = np.linspace(1 - epsilon, 1 + epsilon, 100)

vol = []
potential_atom = []
pressure = []

# Ref is reference, cu is the one strain is applied to
cuRef = bulk("Cu", "fcc", a=3.6, cubic=True)
cuRef.set_calculator(calc)

cu = bulk("Cu", "fcc", a=3.6, cubic=True)
cu.set_calculator(calc)

# For all the strains, apply strain and obtain pot energy & pressure
for i in range(len(strains)):
    cell = cuRef.get_cell()
    cell *= strains[i]
    cu.set_cell(cell, scale_atoms=True)

    # Volume per atom calculation, assuming 4 atoms altogether
    vol.append(((3.6 * strains[i]) ** 3) / 4)

    # Potential energy per atom calculation
    potential = cu.get_potential_energy() / cu.get_global_number_of_atoms()
    potential_atom.append(potential)

    if i == 0:
        min_potential = potential

    if potential < min_potential:
        min_potential = potential
        min_index = i
    # Pressure calculation
    pressure.append(-(1 / 3) * np.trace(cu.get_stress(voigt=False)))

#bulk modulus calculation at eq volume

vol_eq = vol[min_index]
dPdV = (pressure[min_index + 1] - pressure[min_index - 1])/ (vol[min_index + 1] - vol[min_index - 1])

print("Bulk modulus {} GPa".format((- vol_eq * dPdV)/ GPa))

plt.plot(vol, potential_atom)
plt.xlabel('Volume per atom/ Angstroms^3')
plt.ylabel('Potential Energy per atom/ eV')
plt.title("Potential Energy per atom against Volume per atom")
plt.show()

plt.plot(vol, pressure)
plt.xlabel('Volume per atom/ Angstroms^3')
plt.ylabel('Pressure/ eV/ Angstroms^3')
plt.title("Pressure vs Volume per atom")
plt.show()
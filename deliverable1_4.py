import sys
import numpy as np
from ase.build import bulk
from ase.units import GPa

sys.path.append(r'files')
import Morse

calc = Morse.MorsePotential()

shearArray = np.linspace(0.01, 0.05, 50)

shearModArray = []

cuRef = bulk("Cu", "fcc", a = 3.60288, cubic=True)
cuRef.set_calculator(calc)

cu = bulk("Cu", "fcc", a = 3.60288, cubic=True)
cu.set_calculator(calc)

for shear in shearArray:
    cell = cuRef.get_cell()
    cell[0][1] += shear * 3.60288
    cu.set_cell(cell, scale_atoms = True)

    # Shear stress (eV/Å^3), transverse displacement (Å)
    tempStress = cu.get_stress(voigt = False)[0][1]

    # Shear modulus, in GPa
    shearModArray.append((tempStress * 3.6 / (shear * 3.6)) / GPa)

print("Shear modulus: " + str(round(sum(shearModArray) / len(shearModArray), 1)) + "GPa")


# For positive x strain, y and z will be negative
yzStrainArray = np.linspace(0.99, 1, 100)

yzStressArray = []

# 2.617234468937876
cuRef = bulk("Cu", "fcc", a=3.60288, cubic=True)
cuRef.set_calculator(calc)
cu = bulk("Cu", "fcc", a=3.60288, cubic=True)
cu.set_calculator(calc)

for strain in yzStrainArray:
    # Strain x by 1%
    cell = cuRef.get_cell()
    cell[0][0] *= 1.01
    cell[1][1] *= strain
    cell[2][2] *= strain
    cu.set_cell(cell, scale_atoms=True)

    # absolute values, pick one that is the smallest
    yzStressArray.append(abs(cu.get_stress(voigt=False)[1][1]))

# Position in array of lowest stress in yz
minStressIndex = np.argmin(yzStressArray)
poissonRatio = round(100 - (yzStrainArray[minStressIndex] * 100), 3)

print("Poisson ratio: " + str(poissonRatio))
print("Experimental Poisson ratio: 0.355")
absDifference = abs(0.355 - poissonRatio)
print("Percentage difference: " + str(round(absDifference / 0.355 * 100, 1)) + "%")
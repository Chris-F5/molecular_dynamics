import numpy as np

OH_bond_length = 0.1
HH_bond_length = 0.16330

def create_water_molecule(O_pos, charge_orientation, normal):
  np.linalg.norm(charge_orientation)
  HH_axis = np.cross(charge_orientation, normal)
  np.linalg.norm(HH_axis)
  adjacent = (OH_bond_length**2 - (HH_bond_length/2)**2)**(1/2)
  H1_pos = O_pos - charge_orientation * adjacent + HH_axis * HH_bond_length / 2
  H2_pos = O_pos - charge_orientation * adjacent - HH_axis * HH_bond_length / 2
  return [O_pos, H1_pos, H2_pos]

mols = []
mols.append(create_water_molecule([1, 1, 1], np.array([1,0,0]), np.array([0,1,0])))
mols.append(create_water_molecule([1.5, 1, 1], np.array([-1,0,0]), np.array([0,1,0])))

print("SAMPLESYSTEM")
print(f"{len(mols)*3}")
for i,mol in enumerate(mols):
  print(f"{i+1: >5}SOL     OW {i*3+1: >4} {mol[0][0]:>7} {mol[0][1]:>7} {mol[0][2]:>7}")
  print(f"{i+1: >5}SOL    HW1 {i*3+2: >4} {mol[1][0]:>7} {mol[1][1]:>7} {mol[1][2]:>7}")
  print(f"{i+1: >5}SOL    HW2 {i*3+3: >4} {mol[2][0]:>7} {mol[2][1]:>7} {mol[2][2]:>7}")
print("2 2 2") # simulation dims

import sys, math, itertools
import numpy as np
import progressbar

# Model Constants
OH_bond_length = 0.1
HH_bond_length = 0.16330
O_mass = 8
H_mass = 1

simulation_size = 2.3

# Output Files
trajectory_dump_fname = "trajectory.dump"

trajectory_dump_file = open(trajectory_dump_fname, 'w')

# SYSTEM STATE
pos = None
mass = None

def load_system_state(sys_file):
  global pos, last_pos, mass, charge, nmol, force
  simname = sys.stdin.readline().strip()
  print(f"Reading {simname}...", file=sys.stderr)
  natom = int(sys_file.readline())
  assert(natom == 3) # only accept one H2O mol
  pos = np.zeros((natom, 3))
  mass = np.zeros((natom, 3))
  mass[0] = O_mass
  mass[1] = mass[2] = H_mass

  O_mol,O_type,O_n,O_x,O_y,O_z = sys.stdin.readline().split()
  H1_mol,H1_type,H1_n,H1_x,H1_y,H1_z = sys.stdin.readline().split()
  H2_mol,H2_type,H2_n,H2_x,H2_y,H2_z = sys.stdin.readline().split()
  assert(O_mol == H1_mol and O_mol == H2_mol)
  assert(O_type == 'OW' and H1_type == 'HW1' and H2_type == 'HW2')
  pos[0] = np.array([float(c) for c in (O_x, O_y, O_z)])
  pos[1] = np.array([float(c) for c in (H1_x, H1_y, H1_z)])
  pos[2] = np.array([float(c) for c in (H2_x, H2_y, H2_z)])

def dump_timestep(timestep):
  global pos
  print("ITEM: TIMESTEP", file=trajectory_dump_file)
  print(timestep, file=trajectory_dump_file)
  print("ITEM: NUMBER OF ATOMS", file=trajectory_dump_file)
  print(3, file=trajectory_dump_file)
  print("ITEM: BOX BOUNDS pp pp pp", file=trajectory_dump_file)
  print(0, simulation_size, file=trajectory_dump_file)
  print(0, simulation_size, file=trajectory_dump_file)
  print(0, simulation_size, file=trajectory_dump_file)
  print("ITEM: ATOMS id type x y z", file=trajectory_dump_file)
  id = 0
  type = 1
  p = pos[0]
  print(id, type, p[0], p[1], p[2], file=trajectory_dump_file)
  id += 1
  type = 0
  for p in [pos[1], pos[2]]:
    print(id, type, p[0], p[1], p[2], file=trajectory_dump_file)
    id += 1
  trajectory_dump_file.flush()

# SETUP

print("Initialising simulation...", file=sys.stderr)
load_system_state(sys.stdin)

print("Starting simulation...", file=sys.stderr)
dump_timestep(0)
t = 0

# SHAKE

cutoff = 1e-5
bonds = [[0, 1, OH_bond_length], [0,2, OH_bond_length], [1,2,HH_bond_length]]
rs = [pos[i] - pos[j] for i,j,_ in bonds]
iteration = 0
while True:
  iteration += 1
  print(f"iteration: {iteration}")
  worst_constraint = 0
  for (i,j,d),r in zip(bonds, rs):
    r_prime = pos[i] - pos[j]
    rhs = d**2 - np.dot(r_prime, r_prime)
    coefficient = 2 * ((1/mass[i])+(1/mass[j])) * np.dot(r, r_prime)
    g = rhs / coefficient

    pos[i] += g * r / mass[i]
    pos[j] -= g * r / mass[i]

    print(d**2, np.dot(r_prime, r_prime))
    dump_timestep(t)
    t+=1
    relative_error = abs(rhs) / (d**2)
    worst_constraint = max(abs(rhs), relative_error)
  if worst_constraint < cutoff:
    break

trajectory_dump_file.close()

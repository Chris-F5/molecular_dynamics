import numpy as np
import sys

def dump_timestep(timestep, atom_type, atom_position, bounding_box):
  print("ITEM: TIMESTEP")
  print(timestep)
  print("ITEM: NUMBER OF ATOMS")
  print(len(atom_type))
  print("ITEM: BOX BOUNDS pp pp pp")
  print(0, bounding_box[0])
  print(0, bounding_box[1])
  print(0, bounding_box[2])
  print("ITEM: ATOMS id type x y z")
  for id,(type,pos) in enumerate(zip(atom_type, atom_position)):
    print(id, type, pos[0], pos[1], pos[2])

# Units: nm, ns, u

OH_bond_length = 0.1
HH_bond_length = 0.16330
O_mass = 8
H_mass = 1

simulation_size = 2
simulation_duration = 5

# step_size = 5e-7 # realistic size
step_size = 1e-4
num_steps = int(simulation_duration/step_size)

def compute_jacobian(q, out):
  out[0][0] = +(q[0] - q[3])
  out[0][1] = +(q[1] - q[4])
  out[0][2] = +(q[2] - q[5])
  out[0][3] = -(q[0] - q[3])
  out[0][4] = -(q[1] - q[4])
  out[0][5] = -(q[2] - q[5])
  out[0][6] = 0
  out[0][7] = 0
  out[0][8] = 0
  out[1][0] = +(q[0] - q[6])
  out[1][1] = +(q[1] - q[7])
  out[1][2] = +(q[2] - q[8])
  out[1][3] = 0
  out[1][4] = 0
  out[1][5] = 0
  out[1][6] = -(q[0] - q[6])
  out[1][7] = -(q[1] - q[7])
  out[1][8] = -(q[2] - q[8])
  out[2][0] = 0
  out[2][1] = 0
  out[2][2] = 0
  out[2][3] = +(q[3] - q[6])
  out[2][4] = +(q[4] - q[7])
  out[2][5] = +(q[5] - q[8])
  out[2][6] = -(q[3] - q[6])
  out[2][7] = -(q[4] - q[7])
  out[2][8] = -(q[5] - q[8])

mol_mass_vector = np.array([O_mass]*3+[H_mass]*6)
mol_mass_matrix = np.diag(mol_mass_vector)
mol_inverse_mass_matrix = np.linalg.inv(mol_mass_matrix)

# pos, last_pos and force must be vectors of the form (Ox Oy Oz H1x H1y H1z H2x H2y H2z).
def compute_constraint_forces(pos, last_pos, force):
  jacobian = np.zeros((3,9))
  compute_jacobian(pos, jacobian)
  implied_velocity = (pos - last_pos) / step_size
  jacobian_derivative = np.zeros((3,9))
  compute_jacobian(implied_velocity, jacobian_derivative)
  constraint_matrix = - np.matmul(jacobian, np.matmul(mol_inverse_mass_matrix, np.transpose(jacobian)))
  rhs_vector = np.matmul(jacobian_derivative, implied_velocity) + np.matmul(jacobian, np.matmul(mol_inverse_mass_matrix, force))
  constraint_parameters = np.matmul(np.linalg.inv(constraint_matrix), rhs_vector)
  constraint_force = np.matmul(np.transpose(constraint_parameters), jacobian)
  return constraint_force

simname = sys.stdin.readline()
print(f"Reading {simname}", file=sys.stderr)
natoms = int(sys.stdin.readline())
assert(natoms % 3 == 0)
read_pos = []
read_type = []
read_mass = []
for _ in range(natoms//3):
  O_mol,O_type,O_n,O_x,O_y,O_z = sys.stdin.readline().split()
  H1_mol,H1_type,H1_n,H1_x,H1_y,H1_z = sys.stdin.readline().split()
  H2_mol,H2_type,H2_n,H2_x,H2_y,H2_z = sys.stdin.readline().split()
  assert(O_mol == H1_mol and O_mol == H2_mol)
  assert(O_type == 'OW' and H1_type == 'HW1' and H2_type == 'HW2')
  read_pos += [float(c) for c in (O_x,O_y,O_z, H1_x,H1_y,H1_z, H2_x,H2_y,H2_z)]
  read_type += [1,0,0]
  read_mass += list(mol_mass_vector)

atom_type = np.array(read_type)
position = np.array(read_pos)
last_position = np.copy(position)
new_position = np.zeros(position.shape)
mass_vector = np.array(read_mass)

dump_timestep(0, atom_type, np.reshape(position, (-1,3)), [simulation_size for i in range(3)])
for t in range(1, num_steps+1):
  force = np.zeros(position.shape)
  force[3] = 1
  force[9] = 1
  # constraint forces
  for m in range(0, len(position), 9):
    force[m:m+9] += compute_constraint_forces(position[m:m+9], last_position[m:m+9], force[m:m+9])
  # TODO: light restoring force on constraints to account for numerical inaccruacy

  # verlet step
  new_position = 2 * position - last_position + step_size**2 * force / mass_vector
  tmp = last_position
  last_position = position
  position = new_position
  new_position = tmp
  if t % 1000 == 0:
    dump_timestep(t, atom_type, np.reshape(position, (-1,3)), [simulation_size for i in range(3)])

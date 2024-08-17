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

def create_water_molecule(O_pos, charge_orientation, normal):
  np.linalg.norm(charge_orientation)
  HH_axis = np.cross(charge_orientation, normal)
  np.linalg.norm(HH_axis)
  adjacent = (OH_bond_length**2 - (HH_bond_length/2)**2)**(1/2)
  H1_pos = O_pos - charge_orientation * adjacent + HH_axis * HH_bond_length / 2
  H2_pos = O_pos - charge_orientation * adjacent - HH_axis * HH_bond_length / 2
  return [O_pos, H1_pos, H2_pos]

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

O_pos,H1_pos,H2_pos = create_water_molecule([simulation_size/2 for i in range(3)], np.array([1,0,0]), np.array([0,1,1]))
atom_type = np.array([1, 0, 0])
position = np.concatenate((O_pos, H1_pos, H2_pos))
last_position = np.copy(position)
new_position = np.zeros(position.shape)
mass_vector = np.array([O_mass]*3+[H_mass]*6)
mass_matrix = np.diag(mass_vector)
inverse_mass_matrix = np.linalg.inv(mass_matrix)
force = np.zeros(position.shape)
force[0] = 1
implied_velocity = np.zeros(position.shape)
jacobian = np.zeros((3,9))
jacobian_derivative = np.zeros((3,9))

dump_timestep(0, atom_type, np.reshape(position, (3,3)), [simulation_size for i in range(3)])
for t in range(1, num_steps+1):
  # constraint forces
  compute_jacobian(position, jacobian)
  implied_velocity = (position - last_position) / step_size
  compute_jacobian(implied_velocity, jacobian_derivative)
  constraint_matrix = - np.matmul(jacobian, np.matmul(inverse_mass_matrix, np.transpose(jacobian)))
  rhs_vector = np.matmul(jacobian_derivative, implied_velocity) + np.matmul(jacobian, np.matmul(inverse_mass_matrix, force))
  constraint_parameters = np.matmul(np.linalg.inv(constraint_matrix), rhs_vector)
  constraint_force = np.matmul(np.transpose(constraint_parameters), jacobian)
  force += constraint_force

  # verlet step
  new_position = 2 * position - last_position + step_size**2 * force / mass_vector
  last_position = position
  position = new_position
  new_position = np.zeros(position.shape)
  if t % 1000 == 0:
    dump_timestep(t, atom_type, np.reshape(position, (3,3)), [simulation_size for i in range(3)])

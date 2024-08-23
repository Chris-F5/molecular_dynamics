import sys, math, itertools
import numpy as np
import progressbar

# CONSTANTS
# Units: nm, ns, u, e

# Physical Constants
coulomb_constant = 1.38910e8 # non-trivial unit conversion

# Model Constants
OH_bond_length = 0.1
HH_bond_length = 0.16330
O_mass = 8
H_mass = 1
lennard_jones_sigma = 0.316557
lennard_jones_epsilon = 0.650194 * 1e6 # kJ/mol * 1e6 = u(nm)^2/(ns)^2
HW_charge = 0.410
OW_charge = -2 * HW_charge

# Simulation Configuration
simulation_size = 2
simulation_duration = 5e-3
dump_rate = 10

# Numerical Constants
step_size = 5e-7
water_constraint_restore_displacement = 1e10
water_constraint_restore_velocity = 1e10

# Derived Constants
num_steps = int(simulation_duration/step_size)
# Water Constraint Matrices
water_inverse_mass_matrix = np.linalg.inv(np.diag([O_mass]*3 + [H_mass]*6))

# SYSTEM STATE
nmol = 0
pos = None
last_pos = None
O_pos = None
H1_pos = None
H2_pos = None
last_O_pos = None
last_H1_pos = None
last_H2_pos = None
mass = None
charge = None
links = None
def create_atom_views():
  global O_pos, H1_pos, H2_pos, last_O_pos, last_H1_pos, last_H2_pos
  O_pos  = pos[nmol*0:nmol*1]
  H1_pos = pos[nmol*1:nmol*2]
  H2_pos = pos[nmol*2:nmol*3]
  last_O_pos  = last_pos[nmol*0:nmol*1]
  last_H1_pos = last_pos[nmol*1:nmol*2]
  last_H2_pos = last_pos[nmol*2:nmol*3]
def load_system_state(sys_file):
  global pos, last_pos, mass, charge, nmol
  simname = sys.stdin.readline().strip()
  print(f"Reading {simname}...", file=sys.stderr)
  natom = int(sys_file.readline())
  assert(natom % 3 == 0)
  nmol = natom // 3
  pos = np.zeros((natom, 3))
  last_pos = np.zeros((natom, 3))
  mass = np.zeros((natom, 3))
  mass[:nmol] = O_mass
  mass[nmol:] = H_mass
  charge = np.zeros((natom, 3))
  charge[:nmol] = OW_charge
  charge[nmol:] = HW_charge
  create_atom_views()
  for i in range(nmol):
    O_mol,O_type,O_n,O_x,O_y,O_z = sys.stdin.readline().split()
    H1_mol,H1_type,H1_n,H1_x,H1_y,H1_z = sys.stdin.readline().split()
    H2_mol,H2_type,H2_n,H2_x,H2_y,H2_z = sys.stdin.readline().split()
    assert(O_mol == H1_mol and O_mol == H2_mol)
    assert(O_type == 'OW' and H1_type == 'HW1' and H2_type == 'HW2')
    O_pos[i] = np.array([float(c) for c in (O_x, O_y, O_z)])
    H1_pos[i] = np.array([float(c) for c in (H1_x, H1_y, H1_z)])
    H2_pos[i] = np.array([float(c) for c in (H2_x, H2_y, H2_z)])
  last_pos[:] = pos
def dump_timestep(timestep):
  global pos, O_pos, H1_pos, H2_pos
  print("ITEM: TIMESTEP")
  print(timestep)
  print("ITEM: NUMBER OF ATOMS")
  print(nmol*3)
  print("ITEM: BOX BOUNDS pp pp pp")
  print(0, simulation_size)
  print(0, simulation_size)
  print(0, simulation_size)
  print("ITEM: ATOMS id type x y z")
  id = 0
  type = 1
  for p in O_pos:
    print(id, type, p[0], p[1], p[2])
    id += 1
  type = 0
  for p in itertools.chain(H1_pos, H2_pos):
    print(id, type, p[0], p[1], p[2])
    id += 1
def update_links(radius):
  global links
  links = []
  for i in range(nmol):
    O1_pos = O_pos[i]
    for j in range(i+1, nmol):
      best_m = None
      best_r = None
      for mx in [-1,0,1]:
        for my in [-1,0,1]:
          for mz in [-1,0,1]:
            m = np.array([mx, my, mz])
            O2_pos = O_pos[j] + m * simulation_size
            r = np.linalg.norm(O1_pos - O2_pos)
            if best_r == None or r < best_r:
              best_r = r
              best_m = m
      O2_pos = O_pos[j] + best_m * simulation_size
      if best_r <= radius:
        links.append((i, j, best_m))
def enforce_periodic_bc():
  for i in range(nmol):
    for c in range(0,3):
      if O_pos[i][c] < 0:
        O_pos[i][c] += simulation_size
        H1_pos[i][c] += simulation_size
        H2_pos[i][c] += simulation_size
        last_O_pos[i][c] += simulation_size
        last_H1_pos[i][c] += simulation_size
        last_H2_pos[i][c] += simulation_size
      if O_pos[i][c] > simulation_size:
        O_pos[i][c] -= simulation_size
        H1_pos[i][c] -= simulation_size
        H2_pos[i][c] -= simulation_size
        last_O_pos[i][c] -= simulation_size
        last_H1_pos[i][c] -= simulation_size
        last_H2_pos[i][c] -= simulation_size

# FORCES
# Pairwise force functions return foce on first atom in parameters.

def compute_lennard_jones_force(O1_pos, O2_pos):
  r = np.linalg.norm(O2_pos - O1_pos)
  n = 24 * lennard_jones_epsilon * (math.pow(lennard_jones_sigma/r, 6) - math.pow(lennard_jones_sigma/r, 12)) / r
  return (O2_pos - O1_pos) * n / r
def compute_coulomb_force(P1_pos, P2_pos, P1_charge, P2_charge):
  r = np.linalg.norm(P2_pos - P1_pos)
  n = - coulomb_constant * P1_charge * P2_charge / math.pow(r, 2)
  return (P2_pos - P1_pos) * n / r

def compute_jacobian(O_pos, H1_pos, H2_pos):
  jacobian = np.empty((3, 9))
  jacobian[0][0:3] = +(O_pos - H1_pos)
  jacobian[0][3:6] = -(O_pos - H1_pos)
  jacobian[0][6:9] = 0
  jacobian[1][0:3] = +(O_pos - H2_pos)
  jacobian[1][3:6] = 0
  jacobian[1][6:9] = -(O_pos - H2_pos)
  jacobian[2][0:3] = 0
  jacobian[2][3:6] = +(H1_pos - H2_pos)
  jacobian[2][6:9] = -(H1_pos - H2_pos)
  return jacobian
def compute_constraint_function(O_pos, H1_pos, H2_pos):
  c = np.empty(3)
  c[0] = 0.5 * np.linalg.norm(O_pos - H1_pos)
  c[1] = 0.5 * np.linalg.norm(O_pos - H2_pos)
  c[2] = 0.5 * np.linalg.norm(H1_pos - H2_pos)
  return c
def compute_constraint_derivative(jacobian, O_pos, H1_pos, H2_pos, last_O_pos, last_H1_pos, last_H2_pos):
  q = np.concatenate((O_pos, H1_pos, H2_pos))
  q_last = np.concatenate((last_O_pos, last_H1_pos, last_H2_pos))
  vel = (q - q_last)/step_size
  c = np.empty(3)
  return np.matmul(jacobian, vel)
# Returns 3 vectors corresponding to O, H1, H2.
def compute_water_constraint_force(O_force, H1_force, H2_force, O_pos, H1_pos, H2_pos, last_O_pos, last_H1_pos, last_H2_pos):
  O_vel = (O_pos - last_O_pos) / step_size
  H1_vel = (H1_pos - last_H1_pos) / step_size
  H2_vel = (H2_pos - last_H2_pos) / step_size
  fa = np.concatenate((O_force, H1_force, H2_force))
  q = np.concatenate((O_pos, H1_pos, H2_pos))
  qd = np.concatenate((O_vel, H1_vel, H2_vel))
  J = compute_jacobian(O_pos, H1_pos, H2_pos)
  Jd = compute_jacobian(O_vel, H1_vel, H2_vel)
  c = compute_constraint_function(O_pos, H1_pos, H2_pos)
  # cd = compute_constraint_derivative(J, O_pos, H1_pos, H2_pos, last_O_pos, last_H1_pos, last_H2_pos)
  W = water_inverse_mass_matrix
  # ks = water_constraint_restore_displacement
  # kd = water_constraint_restore_velocity
  JT = np.transpose(J)
  lhs_mat = np.matmul(J, np.matmul(W, JT))
  rhs_vec = - np.matmul(J, np.matmul(W, fa)) - np.matmul(Jd, qd) # - ks * c - kd * cd
  #print(- ks * c - kd * cd, rhs_vec, file=sys.stderr)
  constraint_parameters = np.matmul(np.linalg.inv(lhs_mat), rhs_vec)
  constraint_force = np.matmul(np.transpose(constraint_parameters), J)
  constraint_force.resize(3, 3)
  return constraint_force
def constraint_restoring_force(P1, P2, last_P1, last_P2, target_distance):
  displacement = P2 - P1
  distance = np.linalg.norm(displacement)
  direction = displacement / distance
  last_distance = np.linalg.norm(last_P2 - last_P1)

  distance_velocity = (distance - last_distance) / step_size
  return water_constraint_restore_displacement * (distance - target_distance) * direction \
       + water_constraint_restore_velocity * distance_velocity * direction

water_constraint_restore_velocity = 1e5

# SETUP

load_system_state(sys.stdin)
update_links(10)

# TIME PROGRESSION

print("Starting simulation...", file=sys.stderr)
dump_timestep(0)
for t in progressbar.progressbar(range(1, num_steps+1)):
  update_links(10)

  # APPLIED FORCES
  force = np.zeros(pos.shape)
  for i,j,m in links:
    f = compute_lennard_jones_force(pos[i], pos[j] + m * simulation_size)
    force[i] += f
    force[j] -= f
  for i,j,m in links:
    f = compute_coulomb_force(pos[i], pos[j] + m * simulation_size, charge[i], charge[j])
    force[i] += f
    force[j] -= f

  # CONSTRAINT FORCES
  for i in range(nmol):
    O_force = force[i]
    H1_force = force[nmol+i]
    H2_force = force[nmol*2+i]
    O_fc,H1_fc,H2_fc = compute_water_constraint_force(
        O_force, H1_force, H2_force,
        O_pos[i], H1_pos[i], H2_pos[i],
        last_O_pos[i], last_H1_pos[i], last_H2_pos[i],
    )
    force[i] += O_fc
    force[nmol+i] += H1_fc
    force[nmol*2+i] += H2_fc
    O_H1_restoring = constraint_restoring_force(O_pos[i], H1_pos[i], last_O_pos[i], last_H1_pos[i], OH_bond_length)
    O_H2_restoring = constraint_restoring_force(O_pos[i], H2_pos[i], last_O_pos[i], last_H2_pos[i], OH_bond_length)
    H1_H2_restoring = constraint_restoring_force(H1_pos[i], H2_pos[i], last_H1_pos[i], last_H2_pos[i], HH_bond_length)
    force[i] += O_H1_restoring + O_H2_restoring
    force[nmol+i] += H1_H2_restoring - O_H1_restoring
    force[nmol*2+i] += - H1_H2_restoring - O_H2_restoring

  # VERLET INTEGRATION
  new_pos = 2 * pos - last_pos + step_size**2 * force / mass
  last_pos = pos
  pos = new_pos
  create_atom_views()

  enforce_periodic_bc()

  if t % dump_rate == 0:
    dump_timestep(t)

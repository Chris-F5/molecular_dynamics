#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "api.h"
#include "model.h"
#include "utils.h"

const struct particle_params particles_params[] = {
  [PARTICLE_HW] = { 1.0 },
  [PARTICLE_OW] = { 8.0 },
};
const struct non_bonded_interaction_params non_bonded_interactions_params[INTERACTION_COUNT] = {
  [INTERACTION(PARTICLE_HW, PARTICLE_HW)] = { (+0.410)*(+0.410),      0.0,      0.0 },
  [INTERACTION(PARTICLE_OW, PARTICLE_HW)] = { (-0.820)*(+0.410),      0.0,      0.0 },
  [INTERACTION(PARTICLE_OW, PARTICLE_OW)] = { (-0.820)*(-0.820), 0.316557, 0.650194 },
};
const double coulomb_constant = 138.94;
const double non_bonded_cutoff_distance = 1.0;

static void
compute_jacobian(double jacobian[3][9], const double pos[3][3])
{
  displacement(&jacobian[0][0], pos[0], pos[1]);
  displacement(&jacobian[0][3], pos[1], pos[0]);
  jacobian[0][6] = jacobian[0][7] = jacobian[0][8] = 0.0;

  displacement(&jacobian[1][0], pos[0], pos[2]);
  jacobian[1][3] = jacobian[1][4] = jacobian[1][5] = 0.0;
  displacement(&jacobian[1][6], pos[2], pos[0]);

  jacobian[2][0] = jacobian[2][1] = jacobian[2][2] = 0.0;
  displacement(&jacobian[2][3], pos[1], pos[2]);
  displacement(&jacobian[2][6], pos[2], pos[1]);
}

static void
invert_3x3_matrix(double inv[3][3], const double mat[3][3])
{
  double minors[3][3], det;
  int i, j;
  minors[0][0] = mat[1][1] * mat[2][2] - mat[2][1] * mat[1][2];
  minors[0][1] = mat[1][0] * mat[2][2] - mat[2][0] * mat[1][2];
  minors[0][2] = mat[1][0] * mat[2][1] - mat[2][0] * mat[1][1];
  minors[1][0] = mat[0][1] * mat[2][2] - mat[2][1] * mat[0][2];
  minors[1][1] = mat[0][0] * mat[2][2] - mat[2][0] * mat[0][2];
  minors[1][2] = mat[0][0] * mat[2][1] - mat[2][0] * mat[0][1];
  minors[2][0] = mat[0][1] * mat[1][2] - mat[1][1] * mat[0][2];
  minors[2][1] = mat[0][0] * mat[1][2] - mat[1][0] * mat[0][2];
  minors[2][2] = mat[0][0] * mat[1][1] - mat[1][0] * mat[0][1];
  det = mat[0][0] * minors[0][0] + mat[0][1] * minors[0][1] + mat[0][2] * minors[0][2];
  for (i = 0; i < 3; i++) for (j = 0; j < 3; j++)
    inv[i][j] = (det == 0 ? 0 : 1/det) * ((i+j)%2 ? -1 : 1) * minors[j][i];
}

static void
add_spc_restoring_force(const double pos[3][3], const double last_pos[3][3], double force[3][3], double timestep)
{
  int i, j, ax;
  double d[3], target_distance, distance, last_distance, f;
  for (i = 0; i < 3; i++) for (j = i+1; j < 3; j++) {
    target_distance = (i == 0) ? 0.1 : 0.16330; /* TODO: move bond length into params. */
    displacement(d, last_pos[j], pos[i]);
    last_distance = magnitude(d);
    displacement(d, pos[j], pos[i]);
    distance = magnitude(d);
    f = 1e4 * (distance - target_distance) + 1e4 * (distance - last_distance) / timestep;
    for (ax = 0; ax < 3; ax++) {
      force[i][ax] += d[ax] * f / distance;
      force[j][ax] -= d[ax] * f / distance;
    }
  }
}

static void
add_spc_constraint(const double pos[3][3], const double last_pos[3][3], double force[3][3], double timestep)
{
  double vel[3][3], mass;
  double jacobian[3][9], jacobian_derivative[3][9];
  double lhs_matrix[3][3], lhs_inv[3][3];
  double rhs_vector[3], legrangian_multipliers[3];
  int i, j, k;
  for (i = 0; i < 3; i++) for (j = 0; j < 3; j++)
    vel[i][j] = (last_pos[i][j] - pos[i][j]) / timestep;
  compute_jacobian(jacobian, pos);
  compute_jacobian(jacobian_derivative, vel);
  for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) {
    lhs_matrix[i][j] = 0.0;
    for (k = 0; k < 9; k++) {
      mass = particles_params[k/3 ? PARTICLE_OW : PARTICLE_HW].mass;
      lhs_matrix[i][j] += jacobian[i][k] * mass * jacobian[j][k];
    }
  }
  invert_3x3_matrix(lhs_inv, lhs_matrix);
  for (i = 0; i < 3; i++) {
    rhs_vector[i] = 0;
    for (j = 0; j < 9; j++) {
      mass = particles_params[j/3 ? PARTICLE_OW : PARTICLE_HW].mass;
      rhs_vector[i] -= jacobian_derivative[i][j] * vel[j/3][j%3] + jacobian[i][j] * mass * force[j/3][j%3];
    }
  }
  for (i = 0; i < 3; i++) {
    legrangian_multipliers[i] = 0;
    for (j = 0; j < 3; j++)
      legrangian_multipliers[i] += lhs_inv[i][j] * rhs_vector[j];
  }
  for (i = 0; i < 9; i++)
    for (j = 0; j < 3; j++)
      force[i/3][i%3] += legrangian_multipliers[j] * jacobian[j][i];
}

void
enforce_periodic_boundary_conditions(
    double (*p_pos)[3], double (*p_last_pos)[3],
    int water_mol_count, const int (*w_mol)[3],
    double size)
{
  int m, ax;
  for (m = 0; m < water_mol_count; m++) {
    for (ax = 0; ax < 3; ax++) {
      if (p_pos[w_mol[m][0]][ax] > size) {
        p_pos[w_mol[m][0]][ax] -= size;
        p_pos[w_mol[m][1]][ax] -= size;
        p_pos[w_mol[m][2]][ax] -= size;
        if (p_last_pos) {
          p_last_pos[w_mol[m][0]][ax] -= size;
          p_last_pos[w_mol[m][1]][ax] -= size;
          p_last_pos[w_mol[m][2]][ax] -= size;
        }
      } else if (p_pos[w_mol[m][0]][ax] < 0) {
        p_pos[w_mol[m][0]][ax] += size;
        p_pos[w_mol[m][1]][ax] += size;
        p_pos[w_mol[m][2]][ax] += size;
        if (p_last_pos) {
          p_last_pos[w_mol[m][0]][ax] += size;
          p_last_pos[w_mol[m][1]][ax] += size;
          p_last_pos[w_mol[m][2]][ax] += size;
        }
      }
    }
  }
}

/* Assume we are given contiguious molecules of 3 atoms. */
void find_pair_list(
    int particle_count, const double (*p_pos)[3], const int *p_type,
    int interaction_allocated, int *interaction_count, int (*i_particles)[2],
    unsigned char (*i_image), struct non_bonded_interaction_params const **i_params,
    double size)
{
  int pi, pj, ax;
  double d[3];
  *interaction_count = 0;
  for (pi = 0; pi < particle_count; pi++) {
    for (pj = pi-(pi%3)+3; pj < particle_count; pj++) {
      displacement(d, p_pos[pj], p_pos[pi]);
      if (magnitude(d) < non_bonded_cutoff_distance) {
        if (*interaction_count >= interaction_allocated)
          exit(1);
        i_particles[*interaction_count][0] = pi;
        i_particles[*interaction_count][1] = pj;
        i_params[*interaction_count] = &non_bonded_interactions_params[INTERACTION(p_type[pi], p_type[pj])];
        i_image[*interaction_count] = 0;
        for (ax = 0; ax < 3; ax++) {
          if (d[ax]*2 < -size)
            i_image[*interaction_count] |= 0b01 << (ax*2);
          else if (d[ax]*2 > size)
            i_image[*interaction_count] |= 0b10 << (ax*2);
        }
        (*interaction_count)++;
      }
    }
  }
}

double
compute_potential(
    const double (*p_pos)[3],
    int interaction_count, const int (*i_particles)[2], const unsigned char (*i_image),
    const struct non_bonded_interaction_params  *const *i_params,
    double size)
{
  int i;
  double image[3];
  double d[3], r, lj_potential, coulomb_potential, total_potential;
  total_potential = 0;
  for (i = 0; i < interaction_count; i++) {
    image_transform(image, p_pos[i_particles[i][1]], i_image[i], size);
    displacement(d, image, p_pos[i_particles[i][0]]);
    r = magnitude(d);
    lj_potential = 4 * i_params[i]->lennard_jones_epsilon
        * ( pow(i_params[i]->lennard_jones_sigma / r, 12) - pow(i_params[i]->lennard_jones_sigma / r, 6) );
    coulomb_potential = coulomb_constant * i_params[i]->charge / r;
    total_potential += lj_potential + coulomb_potential;
  }
  return total_potential;
}

void add_non_bonded_forces(
    const double (*p_pos)[3], double (*p_force)[3],
    int interaction_count, const int (*i_particles)[2], const unsigned char (*i_image),
    const struct non_bonded_interaction_params  *const *i_params,
    double size)
{
  int i, ax;
  double image[3], d[3];
  double r, lj_magnitude, coulomb_magnitude;
  for (i = 0; i < interaction_count; i++) {
    image_transform(image, p_pos[i_particles[i][1]], i_image[i], size);
    displacement(d, image, p_pos[i_particles[i][0]]);
    r = magnitude(d);
    lj_magnitude = 24 * i_params[i]->lennard_jones_epsilon
        * ( pow(i_params[i]->lennard_jones_sigma / r, 6) - pow(i_params[i]->lennard_jones_sigma / r, 12) )
        / r;
    coulomb_magnitude =  -coulomb_constant * i_params[i]->charge / pow(r, 2);
    for (ax = 0; ax < 3; ax++) {
      p_force[i_particles[i][0]][ax] += d[ax] * (lj_magnitude + coulomb_magnitude) / r;
      p_force[i_particles[i][1]][ax] -= d[ax] * (lj_magnitude + coulomb_magnitude) / r;
    }
  }
}


void add_bonded_forces(
    const double (*p_pos)[3], const double (*p_last_pos)[3], double (*p_force)[3],
    int water_mol_count, const int (*w_mol)[3], double timestep)
{
  int m, atom, ax;
  double pos[3][3];
  double last_pos[3][3];
  double force[3][3];
  for (m = 0; m < water_mol_count; m++) {
    for (atom = 0; atom < 3; atom++) for (ax = 0; ax < 3; ax++) {
      pos[atom][ax] = p_pos[w_mol[m][atom]][ax];
      last_pos[atom][ax] = p_last_pos[w_mol[m][atom]][ax];
      force[atom][ax] = p_force[w_mol[m][atom]][ax];
    }
    add_spc_constraint(pos, last_pos, force, timestep);
    add_spc_restoring_force(pos, last_pos, force, timestep);
    for (atom = 0; atom < 3; atom++) for (ax = 0; ax < 3; ax++)
      p_force[w_mol[m][atom]][ax] = force[atom][ax];
  }
}

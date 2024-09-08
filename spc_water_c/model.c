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
apply_spc_constraint(double pos[3][3], double tolerence)
{
  double c_length[] = {0.1, 0.1, 0.16330}; /* Bond lengths. */
  int c_i[]         = {  0,   0,       1};
  int c_j[]         = {  1,   2,       2};
  double c_displacement[3][3];
  double current_displacement[3];
  int c, ax, iterations;
  double max_err, square_error, err, m, g;
  for (c = 0; c < 3; c++)
    displacement(c_displacement[c], pos[c_i[c]], pos[c_j[c]]);
  iterations = 0;
  do {
    iterations++;
    max_err = 0;
    for (c = 0; c < 3; c++) {
      displacement(current_displacement, pos[c_i[c]], pos[c_j[c]]);
      square_error = pow(c_length[c], 2) - dot(current_displacement, current_displacement);
      m = 1/particles_params[c_i[c] ? PARTICLE_HW : PARTICLE_OW].mass
        + 1/particles_params[PARTICLE_HW].mass;
      g = square_error / (2 * m * dot(c_displacement[c], current_displacement));
      for (ax = 0; ax < 3; ax++) {
        pos[c_i[c]][ax] += g * c_displacement[c][ax] / particles_params[c_i[c] ? PARTICLE_HW : PARTICLE_OW].mass;
        pos[c_j[c]][ax] -= g * c_displacement[c][ax] / particles_params[PARTICLE_HW].mass;
      }
      err = square_error / pow(c_length[c], 2);
      if (err > max_err)
        max_err = err;
    }
  } while(max_err > tolerence);
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


void
apply_bond_constraints(double (*p_pos)[3], int water_mol_count, const int (*w_mol)[3], double tolerence)
{
  int m, atom, ax;
  double pos[3][3];
  for (m = 0; m < water_mol_count; m++) {
    for (atom = 0; atom < 3; atom++) for (ax = 0; ax < 3; ax++)
      pos[atom][ax] = p_pos[w_mol[m][atom]][ax];
    apply_spc_constraint(pos, tolerence);
    for (atom = 0; atom < 3; atom++) for (ax = 0; ax < 3; ax++)
      p_pos[w_mol[m][atom]][ax] = pos[atom][ax];
  }
}

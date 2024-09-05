#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "utils.h"
#include "api.h"
#include "model.h"

#define ALLOCATED_INTERACTIONS 1000000

static double
find_max_force(int particle_count, double (*p_force)[3])
{
  int i;
  double max, force;
  max = 0;
  for (i = 0; i < particle_count; i++) {
    force = magnitude(p_force[i]) ;
    if (force > max)
      max = force;
  }
  return max;
}

int
main(int argc, char **argv)
{
  FILE *traj_out, *stats_out;
  struct gro_structure gro;
  int particle_count;
  double (*p_pos)[3];
  double (*p_force)[3];
  int *p_type;
  int water_mol_count;
  int (*w_mol)[3];
  int interaction_count;
  int (*i_particles)[2];
  unsigned char *i_image;
  struct non_bonded_interaction_params const **i_params;
  double size;

  int t, p, ax;
  double max_force, potential;
  double descent_delta;
  double mass;

  load_gro_structure(stdin, &gro);
  particle_count = gro.particle_count;
  p_pos = xmalloc(particle_count * sizeof(p_pos[0]));
  memcpy(p_pos, gro.p_pos, particle_count * sizeof(p_pos[0]));
  p_type = xmalloc(particle_count * sizeof(p_type[0]));
  memcpy(p_type, gro.p_type, particle_count * sizeof(p_type[0]));
  water_mol_count = gro.water_mol_count;
  w_mol = xmalloc(water_mol_count * sizeof(w_mol[0]));
  memcpy(w_mol, gro.w_mol, water_mol_count * sizeof(w_mol[0]));
  if (gro.dx != gro.dy || gro.dx != gro.dz) {
    fprintf(stderr, "Simulation box must be cube.\n");
    exit(1);
  }
  size = gro.dx;
  free_gro_structure(&gro);

  i_particles = xmalloc(ALLOCATED_INTERACTIONS * sizeof(i_particles[0]));
  i_image = xmalloc(ALLOCATED_INTERACTIONS * sizeof(i_image[0]));
  i_params = xmalloc(ALLOCATED_INTERACTIONS * sizeof(i_params[0]));
  enforce_periodic_boundary_conditions(p_pos, NULL, water_mol_count, w_mol, size);
  find_pair_list(particle_count, p_pos, p_type, ALLOCATED_INTERACTIONS,
      &interaction_count, i_particles, i_image, i_params, size);

  p_force = xmalloc(particle_count * sizeof(p_force[0]));
  memset(p_force, 0, particle_count * sizeof(p_force[0]));
  add_non_bonded_forces(p_pos, p_force, interaction_count, i_particles, i_image, i_params, size);
  add_bonded_forces(p_pos, p_pos, p_force, water_mol_count, w_mol, 1.0);

  descent_delta = 0.01;

  traj_out = fopen("out.traj", "w");
  stats_out = fopen("out.stats", "w");
  fprintf(stats_out, "timestep,potential,max_force\n");
  write_traj_timestep(traj_out, 0, particle_count, p_type, p_pos, size, size, size);
  for (t = 0; t < 800; t++) {
    if (t % 10 == 0) {
      enforce_periodic_boundary_conditions(p_pos, NULL, water_mol_count, w_mol, size);
      find_pair_list(particle_count, p_pos, p_type, ALLOCATED_INTERACTIONS,
          &interaction_count, i_particles, i_image, i_params, size);
    }

    max_force = find_max_force(particle_count, p_force);
    if (max_force > 0) {
      for (p = 0; p < particle_count; p++) for (ax = 0; ax < 3; ax++) {
        mass = particles_params[p_type[p]].mass;
        p_pos[p][ax] = p_pos[p][ax] + descent_delta * (p_force[p][ax] / mass) / max_force;
      }
    }

    memset(p_force, 0, particle_count * sizeof(p_force[0]));
    add_non_bonded_forces(p_pos, p_force, interaction_count, i_particles, i_image, i_params, size);
    add_bonded_forces(p_pos, p_pos, p_force, water_mol_count, w_mol, 1.0);

    if (t % 10 == 0) {
      write_traj_timestep(traj_out, t, particle_count, p_type, p_pos, size, size, size);
      potential = compute_potential(p_pos, interaction_count, i_particles, i_image, i_params, size);
      fprintf(stats_out, "%d,%lf,%lf\n", t, potential, max_force);
    }
  }
  fclose(traj_out);
  fclose(stats_out);

  free(p_pos);
  free(p_type);
  free(w_mol);
  return 0;
}

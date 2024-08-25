#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <regex.h>
#include <errno.h>

#define MAX(a,b) ((a) >= (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define INTERACTION(P1, P2) (MAX(P1, P2) * (MAX(P1, P2) + 1) / 2 + MIN(P1, P2))
#define INTERACTION_COUNT  (PARTICLE_TYPE_COUNT * (PARTICLE_TYPE_COUNT + 1) / 2)

enum ParticleType {
  PARTICLE_HW = 0,
  PARTICLE_OW = 1,
  PARTICLE_TYPE_COUNT = 2,
};

enum integration_type {
  INTEGRATION_VERLET = 0,
  INTEGRATION_STEEPEST_DESCENT = 1,
};

struct particle_params {
  double mass;
} particle_params[] = {
  [PARTICLE_HW] = { 1.0 },
  [PARTICLE_OW] = { 8.0 },
};

struct non_bonded_interaction_params {
  double charge;
  double lennard_jones_sigma;
  double lennard_jones_epsilon;
} non_bonded_interaction_params[INTERACTION_COUNT] = {
  [INTERACTION(PARTICLE_HW, PARTICLE_HW)] = { (+0.410)*(+0.410),      0.0,            0.0 },
  [INTERACTION(PARTICLE_OW, PARTICLE_HW)] = { (-0.820)*(+0.410),      0.0,            0.0 },
  [INTERACTION(PARTICLE_OW, PARTICLE_OW)] = { (-0.820)*(-0.820), 0.316557, 0.650194 * 1e6 },
};

struct simulation_params {
  int integration_type;

  double verlet_step_size;
  double force_cutoff;
  double max_frame_displacement;

  int link_update_frequency;
  int trajectory_dump_frequency;
  int stats_dump_frequency;

  int num_steps;
  double steepest_descent_displacement_end;
};

struct system {
  int timestep;
  double size;

  int particle_count;
  double (*p_pos)[3];
  double (*p_last_pos)[3];
  double (*p_force)[3];
  int *p_type;

  int water_mol_count;
  int (*w_mol)[3];

  int non_bonded_interactions_allocated;
  int non_bonded_interactions_count;
  struct non_bonded_interaction_params *i_params;
  int (*i_particles)[2];

  double steepest_descent_delta;
};

void *
xmalloc(size_t size)
{
  void *p;

  if (!(p = malloc(size))) {
    fprintf(stderr, "malloc failed: %s\n", strerror(errno));
    exit(1);
  }
  return p;
}

void
allocate_system(struct system *sys, int particle_count, int water_mol_count, int non_bonded_interactions_allocated)
{
  sys->particle_count = particle_count;
  sys->p_pos = xmalloc(sys->particle_count * sizeof(sys->p_pos[0]));
  sys->p_last_pos = xmalloc(sys->particle_count * sizeof(sys->p_last_pos[0]));
  sys->p_force = xmalloc(sys->particle_count * sizeof(sys->p_force[0]));
  sys->p_type = xmalloc(sys->particle_count * sizeof(sys->p_type[0]));

  sys->water_mol_count = water_mol_count;
  sys->w_mol = xmalloc(sys->water_mol_count * sizeof(sys->w_mol[0]));

  sys->non_bonded_interactions_allocated = non_bonded_interactions_allocated;
  sys->i_particles = xmalloc(sys->non_bonded_interactions_count * sizeof(sys->i_particles[0]));
}

void
load_gro_structure(FILE *gro_file, struct system *sys)
{
  char line_buffer[1024];
  int p, i, particle_count;
  regex_t r;
  regmatch_t rmatch[7];

  if (fgets(line_buffer, sizeof(line_buffer), gro_file) == NULL)
    goto read_failed;
  line_buffer[strcspn(line_buffer, "\r\n")] = '\0';
  printf("Reading gromacs structure %s...\n", line_buffer);

  if (fgets(line_buffer, sizeof(line_buffer), gro_file) == NULL)
    goto read_failed;
  if (regcomp(&r, "^\\s*([0-9]+)\\s*\r?\n$", REG_EXTENDED))
    goto regex_error;
  if (regexec(&r, line_buffer, 2, rmatch, 0) || rmatch[1].rm_so < 0)
    goto syntax_error;
  regfree(&r);
  line_buffer[rmatch[1].rm_eo] = '\0';
  particle_count = atoi(line_buffer + rmatch[1].rm_so);
  if (particle_count % 3)
    goto unsupported_error;
  allocate_system(sys, particle_count, particle_count/3, 16384);

  if (regcomp(&r, "^\\s*([0-9]+)SOL\\s+(OW|HW1|HW2)\\s+([0-9]+)\\s+([-.0-9]+)\\s+([-.0-9]+)\\s+([-.0-9]+)\\s*\r?\n$", REG_EXTENDED))
    goto regex_error;
  for (p = 0; p < sys->particle_count; p++) {
    if (fgets(line_buffer, sizeof(line_buffer), gro_file) == NULL)
      goto read_failed;
    if (regexec(&r, line_buffer, 7, rmatch, 0) || rmatch[6].rm_so < 0)
      goto syntax_error;
    for (i = 1; i < 7; i++)
      line_buffer[rmatch[i].rm_eo] = '\0';
    if (atoi(line_buffer + rmatch[1].rm_so) != p/3 + 1)
      goto unsupported_error;
    if ( (p % 3 == 0 && strcmp(line_buffer + rmatch[2].rm_so, "OW"))
      || (p % 3 == 1 && strcmp(line_buffer + rmatch[2].rm_so, "HW1"))
      || (p % 3 == 2 && strcmp(line_buffer + rmatch[2].rm_so, "HW2")))
      goto unsupported_error;
    if (atoi(line_buffer + rmatch[3].rm_so) != p + 1)
      goto unsupported_error;
    /* TODO: write particle to system. */
    sys->p_pos[p][0] = sys->p_last_pos[p][0] = atof(line_buffer + rmatch[4].rm_so);
    sys->p_pos[p][1] = sys->p_last_pos[p][1] = atof(line_buffer + rmatch[5].rm_so);
    sys->p_pos[p][2] = sys->p_last_pos[p][2] = atof(line_buffer + rmatch[6].rm_so);
    sys->p_force[p][0] = sys->p_force[p][1] = sys->p_force[p][2] = 0.0;
    sys->p_type[p] = p % 3 ? PARTICLE_HW : PARTICLE_OW;

    if (i % 3 == 0) {
      sys->w_mol[p/3][0] = i;
      sys->w_mol[p/3][1] = i + 1;
      sys->w_mol[p/3][2] = i + 2;
    }
  }

  if (fgets(line_buffer, sizeof(line_buffer), gro_file) == NULL)
    goto read_failed;
  if (regcomp(&r, "^\\s*([.0-9]+)\\s+([.0-9]+)\\s+([.0-9]+)\\s*\r?\n$", REG_EXTENDED))
    goto regex_error;
  if (regexec(&r, line_buffer, 4, rmatch, 0) || rmatch[3].rm_so < 0)
    goto syntax_error;
  regfree(&r);
  for (i = 0; i < 3; i++)
      line_buffer[rmatch[i].rm_eo] = '\0';
  sys->size = atof(line_buffer + rmatch[1].rm_so);
  if (sys->size != atof(line_buffer + rmatch[2].rm_so)
   || sys->size != atof(line_buffer + rmatch[3].rm_so))
    goto unsupported_error;

  sys->timestep = 0;
  sys->steepest_descent_delta = -1.0;
  sys->non_bonded_interactions_count = 0;

  return;
regex_error:
  fflush(stdout);
  fprintf(stderr, "Failed to compile regex for gro structure parsing.\n");
  exit(1);
syntax_error:
  fflush(stdout);
  fprintf(stderr, "Syntax error in gro structure file.\n");
  exit(1);
unsupported_error:
  fflush(stdout);
  fprintf(stderr, "Unsupported gro structure file.\n");
  exit(1);
read_failed:
  fflush(stdout);
  fprintf(stderr, "Error reading from gro structure file.\n");
  exit(1);
}

static void
dump_trajectory_timestep(FILE *traj_file, const struct system *sys)
{
  int i;
  fprintf(traj_file, "ITEM: TIMESTEP\n%d\n", sys->timestep);
  fprintf(traj_file, "ITEM: NUMBER OF ATOMS\n%d\n", sys->particle_count);
  fprintf(traj_file, "ITEM: BOX BOUNDS pp pp pp\n0 %lf\n0 %lf\n0 %lf\n", sys->size, sys->size, sys->size);
  fprintf(traj_file, "ITEM: ATOMS id type x y z\n");
  for (i = -0; i < sys->particle_count; i++)
    fprintf(traj_file, "%d %d %lf %lf %lf\n", i, sys->p_type[i], sys->p_pos[i][0], sys->p_pos[i][1], sys->p_pos[i][2]);
  fflush(traj_file);
}

int
main(int argc, char **argv)
{
  FILE *gro_structure_file;
  struct system sys;

  gro_structure_file = stdin;
  load_gro_structure(gro_structure_file, &sys);
  dump_trajectory_timestep(stdout, &sys);
  return 0;
}

#include <stdio.h>
#include <regex.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "api.h"
#include "utils.h"

void
load_gro_structure(FILE *gro_file, struct gro_structure *gro)
{
  char line_buffer[1024];
  int p, i;
  regex_t r;
  regmatch_t rmatch[7];

  if (fgets(line_buffer, sizeof(line_buffer), gro_file) == NULL)
    goto read_failed;
  line_buffer[strcspn(line_buffer, "\r\n")] = '\0';
  printf("Reading gromacs structure %s...\n", line_buffer);
  strncpy(gro->name, line_buffer, sizeof(gro->name) - 1);
  gro->name[sizeof(gro->name)-1] = '\0';

  if (fgets(line_buffer, sizeof(line_buffer), gro_file) == NULL)
    goto read_failed;
  if (regcomp(&r, "^\\s*([0-9]+)\\s*\r?\n$", REG_EXTENDED))
    goto regex_error;
  if (regexec(&r, line_buffer, 2, rmatch, 0) || rmatch[1].rm_so < 0)
    goto syntax_error;
  regfree(&r);
  line_buffer[rmatch[1].rm_eo] = '\0';
  gro->particle_count = atoi(line_buffer + rmatch[1].rm_so);
  if (gro->particle_count % 3)
    goto unsupported_error;

  gro->p_pos = xmalloc(gro->particle_count * sizeof(gro->p_pos[0]));
  gro->p_type = xmalloc(gro->particle_count * sizeof(gro->p_type[0]));

  gro->water_mol_count = gro->particle_count / 3;;
  gro->w_mol = xmalloc(gro->water_mol_count * sizeof(gro->w_mol[0]));

  if (regcomp(&r, "^\\s*([0-9]+)SOL\\s+(OW|HW1|HW2)\\s+([0-9]+)\\s+([-.0-9]+)\\s+([-.0-9]+)\\s+([-.0-9]+)\\s*\r?\n$", REG_EXTENDED))
    goto regex_error;
  for (p = 0; p < gro->particle_count; p++) {
    if (fgets(line_buffer, sizeof(line_buffer), gro_file) == NULL)
      goto read_failed;
    if (regexec(&r, line_buffer, 7, rmatch, 0) || rmatch[6].rm_so < 0)
      goto syntax_error;
    for (i = 1; i <= 6; i++)
      line_buffer[rmatch[i].rm_eo] = '\0';
    if (atoi(line_buffer + rmatch[1].rm_so) != p/3 + 1)
      goto unsupported_error;
    if ( (p % 3 == 0 && strcmp(line_buffer + rmatch[2].rm_so, "OW"))
      || (p % 3 == 1 && strcmp(line_buffer + rmatch[2].rm_so, "HW1"))
      || (p % 3 == 2 && strcmp(line_buffer + rmatch[2].rm_so, "HW2")))
      goto unsupported_error;
    if (atoi(line_buffer + rmatch[3].rm_so) != p + 1)
      goto unsupported_error;
    /* TODO: write particle to grotem. */
    gro->p_pos[p][0] = atof(line_buffer + rmatch[4].rm_so);
    gro->p_pos[p][1] = atof(line_buffer + rmatch[5].rm_so);
    gro->p_pos[p][2] = atof(line_buffer + rmatch[6].rm_so);
    gro->p_type[p] = p % 3 ? PARTICLE_HW : PARTICLE_OW;

    if (p % 3 == 0) {
      gro->w_mol[p/3][0] = p;
      gro->w_mol[p/3][1] = p + 1;
      gro->w_mol[p/3][2] = p + 2;
    }
  }

  if (fgets(line_buffer, sizeof(line_buffer), gro_file) == NULL)
    goto read_failed;
  if (regcomp(&r, "^\\s*([.0-9]+)\\s+([.0-9]+)\\s+([.0-9]+)\\s*\r?\n$", REG_EXTENDED))
    goto regex_error;
  if (regexec(&r, line_buffer, 4, rmatch, 0) || rmatch[3].rm_so < 0)
    goto syntax_error;
  regfree(&r);
  for (i = 1; i <= 3; i++)
      line_buffer[rmatch[i].rm_eo] = '\0';
  gro->dx = atof(line_buffer + rmatch[1].rm_so);
  gro->dy = atof(line_buffer + rmatch[2].rm_so);
  gro->dz = atof(line_buffer + rmatch[3].rm_so);

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

void
free_gro_structure(struct gro_structure *gro)
{
  free(gro->p_type);
  free(gro->p_pos);
  free(gro->w_mol);
}

void write_traj_timestep(FILE *file, int timestep,
    int particle_count, const int *p_type, const double (*p_pos)[3],
    double dx, double dy, double dz)
{
  int i;
  fprintf(file, "ITEM: TIMESTEP\n%d\n", timestep);
  fprintf(file, "ITEM: NUMBER OF ATOMS\n%d\n", particle_count);
  fprintf(file, "ITEM: BOX BOUNDS pp pp pp\n0 %lf\n0 %lf\n0 %lf\n", dx, dy, dz);
  fprintf(file, "ITEM: ATOMS id type x y z\n");
  for (i = -0; i < particle_count; i++)
    fprintf(file, "%d %d %lf %lf %lf\n", i, p_type[i], p_pos[i][0], p_pos[i][1], p_pos[i][2]);
  fflush(file);
}

#include <stdio.h>
#include <math.h>
#include <string.h>

#define MAX_PARTICLE_COUNT 1000
#define DALTON_MASS 1.66053906892e-27
#define BOLTZMANN_CONSTANT 1.380649e-23
#define LENNARD_JONES_EPSILON (10.22*BOLTZMANN_CONSTANT)
#define LENNARD_JONES_SIGMA (2.556e-10)

#define SIMULATION_DURATION 1e-11
const double tick_duration = 1e-14;

struct vec2 {
  double x, y;
};

struct system {
  long tick;
  long p_count;
  struct vec2 p_pos[MAX_PARTICLE_COUNT];
  struct vec2 p_last_pos[MAX_PARTICLE_COUNT];
  double p_mass[MAX_PARTICLE_COUNT];
  struct vec2 p_force[MAX_PARTICLE_COUNT];
};

void
compute_force(struct system *sys)
{
  /* Lennard-Jones Potential Force */
  long i, j;
  double magnitude, dx, dy, r;
  for (i = 0; i < sys->p_count; i++) {
    for (j = i + 1; j < sys->p_count; j++) {
      dx = sys->p_pos[i].x - sys->p_pos[j].x;
      dy = sys->p_pos[i].y - sys->p_pos[j].y;
      r = sqrt(pow(dx, 2) + pow(dy, 2));
      magnitude
          = 24 * LENNARD_JONES_EPSILON
          * (2*pow(LENNARD_JONES_SIGMA/r, 12) - pow(LENNARD_JONES_SIGMA/r, 6)) / r;
      sys->p_force[i].x = magnitude * dx / r;
      sys->p_force[i].y = magnitude * dy / r;
      sys->p_force[j].x = - magnitude * dx / r;
      sys->p_force[j].y = - magnitude * dy / r;
    }
  }
}

double
compute_system_kenetic_energy(struct system *sys)
{
  long i;
  double speed_squared, total;
  total = 0;
  for (i = 0; i < sys->p_count; i++) {
    speed_squared = pow((sys->p_last_pos[i].x - sys->p_pos[i].x)/tick_duration, 2)
                  + pow((sys->p_last_pos[i].y - sys->p_pos[i].y)/tick_duration, 2);
    total += 0.5 * sys->p_mass[i] * speed_squared;
  }
  return total;
}

double
compute_system_potential_energy(struct system *sys)
{
  long i, j;
  double dx, dy, r, total;
  total = 0;
  for (i = 0; i < sys->p_count; i++) {
    for (j = i + 1; j < sys->p_count; j++) {
      dx = sys->p_pos[i].x - sys->p_pos[j].x;
      dy = sys->p_pos[i].y - sys->p_pos[j].y;
      r = sqrt(pow(dx, 2) + pow(dy, 2));
      total += 4 * LENNARD_JONES_EPSILON * (pow(LENNARD_JONES_SIGMA/r, 12) - pow(LENNARD_JONES_SIGMA/r, 6))
             + LENNARD_JONES_EPSILON;
    }
  }
  return total;
}

void
verlet_tick(struct system *sys)
{
  double t;
  for (long i = 0; i < sys->p_count; i++) {
    t = sys->p_pos[i].x;
    //printf("%e * %e * %e / %e = %e\n", tick_duration, tick_duration, sys->p_force[i].x, sys->p_mass[i], tick_duration * tick_duration * sys->p_force[i].x / sys->p_mass[i]);
    sys->p_pos[i].x
        = 2 * sys->p_pos[i].x - sys->p_last_pos[i].x
        + tick_duration * tick_duration * sys->p_force[i].x / sys->p_mass[i];
    sys->p_last_pos[i].x = t;
    t = sys->p_pos[i].y;
    sys->p_pos[i].y
        = 2 * sys->p_pos[i].y - sys->p_last_pos[i].y
        + tick_duration * tick_duration * sys->p_force[i].y / sys->p_mass[i];
    sys->p_last_pos[i].y = t;
  }
}

void
log_events(struct system *sys)
{
  long i;
  double ke, pe;
  ke = compute_system_kenetic_energy(sys);
  pe = compute_system_potential_energy(sys);
  printf("{\"type\":\"energy\",\"time\":%e,\"kenetic\":%e,\"potential\":%e}\n",
      sys->tick * tick_duration, ke, pe);
  for (i = 0; i < sys->p_count; i++) {
    printf("{\"type\":\"pos\", \"id\":%ld, \"time\":%e, \"x\":%e, \"y\":%e}\n",
        i, sys->tick * tick_duration, sys->p_pos[i].x, sys->p_pos[i].y);
    printf("{\"type\":\"force\", \"id\":%ld, \"time\":%e, \"x\":%e, \"y\":%e}\n",
        i, sys->tick * tick_duration, sys->p_force[i].x, sys->p_force[i].y);
  }
}

void
tick_system(struct system *sys)
{
  verlet_tick(sys);
  compute_force(sys);
  log_events(sys);
  sys->tick++;
}

int
main(int argc, char **argv)
{
  struct system sys;
  sys.tick = 0;
  sys.p_count = 2;
  sys.p_pos[0].x = sys.p_last_pos[0].x = 0;
  sys.p_pos[0].y = sys.p_last_pos[0].y = 0;
  sys.p_pos[1].x = sys.p_last_pos[1].x = 4e-10;
  sys.p_pos[1].y = sys.p_last_pos[1].y = 0;
  sys.p_mass[0] = 4 * DALTON_MASS;
  sys.p_mass[1] = 4 * DALTON_MASS;
  memset(&sys.p_force, 0, sizeof(sys.p_force));
  for (long t = 0; t * tick_duration < SIMULATION_DURATION; t++)
    tick_system(&sys);
  return 0;
}

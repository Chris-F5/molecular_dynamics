#include <stdio.h>
#include <math.h>
#include <string.h>

#define MAX_CLUSTER_OCCUPANCY 1000

#define DALTON_MASS 1.66053906892e-27
#define BOLTZMANN_CONSTANT 1.380649e-23
#define LENNARD_JONES_EPSILON (10.22*BOLTZMANN_CONSTANT)
#define LENNARD_JONES_SIGMA (2.556e-10)

#define SIMULATION_DURATION 1e-11
#define TICK_DURATION 1e-14

#define SIMULATION_SIZE 1e-8

struct vec2 {
  double x, y;
};

struct cluster {
  long particle_count;
  long p_id[MAX_CLUSTER_OCCUPANCY];
  struct vec2 p_pos[MAX_CLUSTER_OCCUPANCY];
  struct vec2 p_lpos[MAX_CLUSTER_OCCUPANCY];
  struct vec2 p_force[MAX_CLUSTER_OCCUPANCY];
  double p_mass[MAX_CLUSTER_OCCUPANCY];
};

struct system {
  long tick;
  struct cluster cluster;

  int link_count;
  struct vec2 link_offset[4];
};

void
lennard_jones_force(double dx, double dy, struct vec2 *a_force, struct vec2 *b_force)
{
  double r, magnitude;
  r = sqrt(pow(dx, 2) + pow(dy, 2));
  magnitude
      = 24 * LENNARD_JONES_EPSILON
      * (2*pow(LENNARD_JONES_SIGMA/r, 12) - pow(LENNARD_JONES_SIGMA/r, 6)) / r;
  a_force->x -= magnitude * dx / r;
  a_force->y -= magnitude * dy / r;
  b_force->x += magnitude * dx / r;
  b_force->y += magnitude * dy / r;
}

double
lennard_jones_potential(double dx, double dy)
{
  double r;
  r = sqrt(pow(dx, 2) + pow(dy, 2));
  return 4 * LENNARD_JONES_EPSILON * (pow(LENNARD_JONES_SIGMA/r, 12) - pow(LENNARD_JONES_SIGMA/r, 6))
       + LENNARD_JONES_EPSILON;
}

void
intracluster_force(struct cluster *cluster)
{
  long i, j;
  double dx, dy, magnitude;
  for (i = 0; i < cluster->particle_count; i++)
    cluster->p_force[i].x = cluster->p_force[i].y = 0.0;
  for (i = 0; i < cluster->particle_count; i++) {
    for (j = i + 1; j < cluster->particle_count; j++) {
      dx = cluster->p_pos[j].x - cluster->p_pos[i].x;
      dy = cluster->p_pos[j].y - cluster->p_pos[i].y;
      lennard_jones_force(dx, dy, &cluster->p_force[i], &cluster->p_force[j]);
    }
  }
}

void
intercluster_force(struct cluster *base, struct cluster *other, struct vec2 offset)
{
  long i, j;
  double dx, dy;
  for (i = 0; i < base->particle_count; i++) {
    for (j = 0; j < other->particle_count; j++) {
      dx = offset.x + other->p_pos[j].x - base->p_pos[i].x;
      dy = offset.y + other->p_pos[j].y - base->p_pos[i].y;
      lennard_jones_force(dx, dy, &base->p_force[i], &other->p_force[j]);
    }
  }
}

double
intracluster_potential(const struct cluster *cluster)
{
  long i, j;
  double dx, dy, magnitude, total;
  total = 0.0;
  for (i = 0; i < cluster->particle_count; i++) {
    for (j = i + 1; j < cluster->particle_count; j++) {
      dx = cluster->p_pos[j].x - cluster->p_pos[i].x;
      dy = cluster->p_pos[j].y - cluster->p_pos[i].y;
      total += lennard_jones_potential(dx ,dy);
    }
  }
  return total;
}

double
intercluster_potential(const struct cluster *base, const struct cluster *other, struct vec2 offset)
{
  long i, j;
  double dx, dy, total;
  total = 0.0;
  for (i = 0; i < base->particle_count; i++) {
    for (j = 0; j < other->particle_count; j++) {
      dx = offset.x + other->p_pos[j].x - base->p_pos[i].x;
      dy = offset.y + other->p_pos[j].y - base->p_pos[i].y;
      total += lennard_jones_potential(dx, dy);
    }
  }
  return total;
}

double
cluster_kenetic_energy(const struct cluster *cluster)
{
  long i;
  double speed_squared, total;
  total = 0.0;
  for (i = 0; i < cluster->particle_count; i++) {
    speed_squared = pow((cluster->p_lpos[i].x - cluster->p_pos[i].x)/TICK_DURATION, 2)
                  + pow((cluster->p_lpos[i].y - cluster->p_pos[i].y)/TICK_DURATION, 2);
    total += 0.5 * cluster->p_mass[i] * speed_squared;
  }
  return total;
}

double
system_potential_energy(const struct system *sys)
{
  double total;
  int l;
  total = 0.0;
  total += intracluster_potential(&sys->cluster);
  for (l = 0; l < sys->link_count; l++)
    total += intercluster_potential(&sys->cluster, &sys->cluster, sys->link_offset[l]);
  return total;
}

void
verlet_tick(struct cluster *cluster)
{
  long i;
  double t;
  for (i = 0; i < cluster->particle_count; i++) {
    t = cluster->p_pos[i].x;
    cluster->p_pos[i].x
        = 2 * cluster->p_pos[i].x - cluster->p_lpos[i].x
        + TICK_DURATION * TICK_DURATION * cluster->p_force[i].x / cluster->p_mass[i];
    cluster->p_lpos[i].x = t;
    t = cluster->p_pos[i].y;
    cluster->p_pos[i].y
        = 2 * cluster->p_pos[i].y - cluster->p_lpos[i].y
        + TICK_DURATION * TICK_DURATION * cluster->p_force[i].y / cluster->p_mass[i];
    cluster->p_lpos[i].y = t;
  }
}

void
log_events(struct system *sys)
{
  int c;
  long i;
  double ke, pe;
  ke = cluster_kenetic_energy(&sys->cluster);
  pe = system_potential_energy(sys);
  printf("{\"type\":\"energy\",\"time\":%e,\"kenetic\":%e,\"potential\":%e}\n",
      sys->tick * TICK_DURATION, ke, pe);
  for (i = 0; i < sys->cluster.particle_count; i++) {
    printf("{\"type\":\"pos\", \"id\":%ld, \"time\":%e, \"x\":%e, \"y\":%e}\n",
        sys->cluster.p_id[i], sys->tick * TICK_DURATION, sys->cluster.p_pos[i].x, sys->cluster.p_pos[i].y);
    printf("{\"type\":\"force\", \"id\":%ld, \"time\":%e, \"x\":%e, \"y\":%e}\n",
        sys->cluster.p_id[i], sys->tick * TICK_DURATION, sys->cluster.p_force[i].x, sys->cluster.p_force[i].y);
  }
}

void
tick_system(struct system *sys)
{
  int l;
  long i;
  verlet_tick(&sys->cluster);
  for (i = 0; i < sys->cluster.particle_count; i++) {
    if (sys->cluster.p_pos[i].x > SIMULATION_SIZE) {
      sys->cluster.p_pos[i].x -= SIMULATION_SIZE;
      sys->cluster.p_lpos[i].x -= SIMULATION_SIZE;
    } else if (sys->cluster.p_pos[i].x < 0) {
      sys->cluster.p_pos[i].x += SIMULATION_SIZE;
      sys->cluster.p_lpos[i].x += SIMULATION_SIZE;
    }
    if (sys->cluster.p_pos[i].y > SIMULATION_SIZE) {
      sys->cluster.p_pos[i].y -= SIMULATION_SIZE;
      sys->cluster.p_lpos[i].y -= SIMULATION_SIZE;
    } else if (sys->cluster.p_pos[i].y < 0) {
      sys->cluster.p_pos[i].y += SIMULATION_SIZE;
      sys->cluster.p_lpos[i].y += SIMULATION_SIZE;
    }
  }
  intracluster_force(&sys->cluster);
  for (l = 0; l < sys->link_count; l++)
    intercluster_force(&sys->cluster, &sys->cluster, sys->link_offset[l]);
  log_events(sys);
  sys->tick++;
}

int
main(int argc, char **argv)
{
  struct system sys;
  sys.tick = 0;
  sys.link_count = 4;
  sys.link_offset[0].x = SIMULATION_SIZE;
  sys.link_offset[0].y = 0.0;
  sys.link_offset[1].x = 0.0;
  sys.link_offset[1].y = SIMULATION_SIZE;
  sys.link_offset[2].x = SIMULATION_SIZE;
  sys.link_offset[2].y = SIMULATION_SIZE;
  sys.link_offset[3].x = SIMULATION_SIZE;
  sys.link_offset[3].y = -SIMULATION_SIZE;
  sys.cluster.particle_count = 2;
  sys.cluster.p_id[0] = 0;
  sys.cluster.p_id[1] = 1;
  sys.cluster.p_pos[0].x = sys.cluster.p_lpos[0].x = 1e-10;
  sys.cluster.p_pos[0].y = sys.cluster.p_lpos[0].y = 1e-10;
  sys.cluster.p_pos[1].x = sys.cluster.p_lpos[1].x = 5e-10;
  sys.cluster.p_pos[1].y = sys.cluster.p_lpos[1].y = 1e-10;
  sys.cluster.p_lpos[0].x -= 2e3 * TICK_DURATION;
  sys.cluster.p_lpos[1].x -= 2e3 * TICK_DURATION;
  sys.cluster.p_mass[0] = 4 * DALTON_MASS;
  sys.cluster.p_mass[1] = 4 * DALTON_MASS;
  memset(&sys.cluster.p_force, 0, sizeof(sys.cluster.p_force));
  for (long t = 0; t * TICK_DURATION < SIMULATION_DURATION; t++)
    tick_system(&sys);
  return 0;
}

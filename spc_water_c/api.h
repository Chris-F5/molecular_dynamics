enum ParticleType {
  PARTICLE_HW = 0,
  PARTICLE_OW = 1,
  PARTICLE_TYPE_COUNT = 2,
};

struct gro_structure {
  char name[256];

  int particle_count;
  int *p_type;
  double (*p_pos)[3];

  int water_mol_count;
  int (*w_mol)[3];

  double dx, dy, dz;
};

void load_gro_structure(FILE *gro_file, struct gro_structure *gro);
void free_gro_structure(struct gro_structure *gro);
void write_traj_timestep(FILE *file, int timestep,
    int particle_count, const int *p_type, const double (*p_pos)[3],
    double dx, double dy, double dz);

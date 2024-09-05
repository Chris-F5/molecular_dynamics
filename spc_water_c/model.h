#define MAX(a,b) ((a) >= (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define INTERACTION(P1, P2) (MAX(P1, P2) * (MAX(P1, P2) + 1) / 2 + MIN(P1, P2))
#define INTERACTION_COUNT  (PARTICLE_TYPE_COUNT * (PARTICLE_TYPE_COUNT + 1) / 2)

struct particle_params {
  double mass;
};
struct non_bonded_interaction_params {
  double charge;
  double lennard_jones_sigma;
  double lennard_jones_epsilon;
};
extern const struct particle_params particles_params[];
extern const struct non_bonded_interaction_params non_bonded_interactions_params[];
extern const double coulomb_constant;
extern const double non_bonded_cutoff_distance;

void enforce_periodic_boundary_conditions(
    double (*p_pos)[3], double (*p_last_pos)[3],
    int water_mol_count, const int (*w_mol)[3],
    double size);
void find_pair_list(
    int particle_count, const double (*p_pos)[3], const int *p_type,
    int interaction_allocated, int *interaction_count, int (*i_particles)[2],
    unsigned char (*i_image), struct non_bonded_interaction_params const **i_params,
    double size);
double compute_potential(
    const double (*p_pos)[3],
    int interaction_count, const int (*i_particles)[2], const unsigned char (*i_image),
    const struct non_bonded_interaction_params  *const *i_params,
    double size);
void add_non_bonded_forces(
    const double (*p_pos)[3], double (*p_force)[3],
    int interaction_count, const int (*i_particles)[2], const unsigned char (*i_image),
    const struct non_bonded_interaction_params  *const *i_params,
    double size);
void add_bonded_forces(
    const double (*p_pos)[3], const double (*p_last_pos)[3], double (*p_force)[3],
    int water_mol_count, const int (*w_mol)[3], double timestep);

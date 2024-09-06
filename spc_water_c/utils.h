void *xmalloc(size_t size);

static inline double
distance(const double a_pos[3], const double b_pos[3])
{
  return sqrt( pow(a_pos[0] - b_pos[0], 2)
             + pow(a_pos[1] - b_pos[1], 2)
             + pow(a_pos[2] - b_pos[2], 2) );
}
static inline double
magnitude(const double vec[3])
{
  return sqrt(pow(vec[0], 2) + pow(vec[1], 2) + pow(vec[2], 2) );
}
static inline void
displacement(double out[3], const double dst_pos[3], const double src_pos[3])
{
  out[0] = dst_pos[0] - src_pos[0];
  out[1] = dst_pos[1] - src_pos[1];
  out[2] = dst_pos[2] - src_pos[2];
}
static inline void
image_transform(double out[3], const double pos[3], unsigned char image, double size)
{
  int ax;
  for (ax = 0; ax < 3; ax++) {
    out[ax] = pos[ax];
    if ( (image >> (ax*2)) & 0b01 )
      out[ax] += size;
    if ( (image >> (ax*2)) & 0b10 )
      out[ax] -= size;
  }
}

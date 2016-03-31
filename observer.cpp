#include "observer.hpp"
#include "parameter.hpp"
#include <numeric>

std::string Observer::fid_to_name(const int fid) const {
  switch (fid) {
  case RadGyr:
    return "./rg.txt";
  case SdRad:
    return "./sdrad.txt";
  case Tempera:
    return "./tempera.txt";
  case Config:
    return "./ptcl_config.xyz";
  default:
    std::cerr << "Do not know fid "  << fid << std::endl;
    std::exit(1);
  }
}

void Observer::calculate_rg(const double3* pos,
                            const int time) {
  const double3 e2e = {pos[0].x - pos[Parameter::MOL_NUM - 1].x,
                       pos[0].y - pos[Parameter::MOL_NUM - 1].y,
                       pos[0].z - pos[Parameter::MOL_NUM - 1].z};
  const double Rg2 = e2e.x * e2e.x + e2e.y * e2e.y + e2e.z * e2e.z;
  const double Rg  = std::sqrt(Rg2);
  fprintf(fp[RadGyr], "%d %g \n", time, Rg);
}

void Observer::calculate_dev_rad(const double3* pos,
                                 const int time) {
  double3 cm_pos = {0.0, 0.0, 0.0};
  double inv_molnum = 1.0 / static_cast<double>(Parameter::MOL_NUM);
  for (int i = 0; i < Parameter::MOL_NUM; i++) {
    cm_pos.x += pos[i].x;
    cm_pos.y += pos[i].y;
    cm_pos.z += pos[i].z;
  }
  cm_pos.x *= inv_molnum;
  cm_pos.y *= inv_molnum;
  cm_pos.z *= inv_molnum;
  
  double sd_rad = 0.0;
  for (int i = 0; i < Parameter::MOL_NUM; i++) {
    const double3 cm2pos = {pos[i].x - cm_pos.x,
                            pos[i].y - cm_pos.y,
                            pos[i].z - cm_pos.z};
    sd_rad += cm2pos.x * cm2pos.x + cm2pos.y * cm2pos.y + cm2pos.z * cm2pos.z;
  }
  sd_rad *= inv_molnum;
  fprintf(fp[SdRad], "%d %g \n", time, sd_rad);
}

void Observer::calculate_tempera(const double3* vel) {
  double sum_vel = std::accumulate(vel, vel + Parameter::MOL_NUM, 0.0,
                                   [](const double sum, const double3& elem){
                                     return sum + elem * elem;
                                   });
  sum_vel /= 3.0 * Parameter::MOL_NUM;
  fprintf(fp[Tempera], "%g \n", sum_vel);
}

void Observer::print_config(const double3* pos, const int time) {
  fprintf(fp[Config] , "%d\n", Parameter::MOL_NUM);
  fprintf(fp[Config], "time %d\n", time);
  for (int i = 0; i < Parameter::MOL_NUM; i++)
    fprintf(fp[Config], "O %.10g %.10g %.10g\n", pos[i].x, pos[i].y, pos[i].z);
}

#include "f_calculator.hpp"
#include "mknearlist.hpp"
#include <iostream>
#include <cstdlib>
#include <limits>

void initialize(double3* pos,
		double3* vel,
		double3* force,
		Plist& pair_list)
{
  const double3 vec3NaN = {std::numeric_limits<double>::quiet_NaN(),
			   std::numeric_limits<double>::quiet_NaN(),
			   std::numeric_limits<double>::quiet_NaN()};
  for(int i = 0; i < Parameter::MOL_NUM; i++){
    pos[i] = vel[i] = vec3NaN;
    force[i].x = force[i].y = force[i].z = 0.0;
  }
  
  const int est_mol_num = static_cast<int>(Parameter::MOL_NUM * 0.3);
  for(int i = 0; i < Parameter::MOL_NUM; i++){
    pair_list[i].resize(est_mol_num, 0);
  }
}

void create_init_config(double3* pos, 
			double3* vel,
			PRNG& prng)
{
  const double T = Parameter::tempera;
  for(int i = 0; i < Parameter::MOL_NUM; i++){
    pos[i].x = Parameter::bond * (i + 0.5); pos[i].y = pos[i].z = 0.0;
    vel[i].x = prng.Normal(0, 0.0, std::sqrt(T)); vel[i].y = prng.Normal(0, 0.0, std::sqrt(T)); vel[i].z = prng.Normal(0, 0.0, std::sqrt(T));
  }
  
  double init_tempera = 0.0;
  for(int i = 0; i < Parameter::MOL_NUM; i++){
    init_tempera += vel[i].x * vel[i].x + vel[i].y * vel[i].y + vel[i].z * vel[i].z;
  }
  init_tempera /= (3.0 * Parameter::MOL_NUM);
  printf("initial temperature is %g \n", init_tempera);
}

void clear_force(double3* force)
{
  for(int i = 0; i < Parameter::MOL_NUM; i++){
    force[i].x = force[i].y = force[i].z = 0.0;
  }
}

void update_pos_vel(double3* pos,
		    double3* vel,
		    const double3* force)
{
  const double cf0 = (1.0 - 0.5 * Parameter::gamm * Parameter::dt);
  for(int i = 0; i < Parameter::MOL_NUM; i++){
    //vel(t) -> vel(t + dt/2)
    vel[i].x = cf0 * vel[i].x + Parameter::half_dt_mass * force[i].x;
    vel[i].y = cf0 * vel[i].y + Parameter::half_dt_mass * force[i].y;
    vel[i].z = cf0 * vel[i].z + Parameter::half_dt_mass * force[i].z;

    //pos(t) -> pos(t + dt)
    pos[i].x += vel[i].x * Parameter::dt;
    pos[i].y += vel[i].y * Parameter::dt;
    pos[i].z += vel[i].z * Parameter::dt;
  }
}

void update_vel(double3* vel,
		const double3* force)
{
  for(int i = 0; i < Parameter::MOL_NUM; i++){
    vel[i].x = Parameter::time_cf * (vel[i].x + Parameter::half_dt_mass * force[i].x);
    vel[i].y = Parameter::time_cf * (vel[i].y + Parameter::half_dt_mass * force[i].y);
    vel[i].z = Parameter::time_cf * (vel[i].z + Parameter::half_dt_mass * force[i].z);
  }
}

void calculate_rg(const double3* pos,
		  const int time,
		  FILE* fp)
{
  const double3 e2e = {pos[0].x - pos[Parameter::MOL_NUM - 1].x,
		       pos[0].y - pos[Parameter::MOL_NUM - 1].y,
		       pos[0].z - pos[Parameter::MOL_NUM - 1].z};
  const double Rg2 = e2e.x * e2e.x + e2e.y * e2e.y + e2e.z * e2e.z;
  const double Rg  = std::sqrt(Rg2);
  fprintf(fp, "%d %g \n", time, Rg);
}

void calculate_dev_rad(const double3* pos,
		       const int time,
		       FILE* fp)
{
  double3 cm_pos = {0.0, 0.0, 0.0};
  double inv_molnum = 1.0 / static_cast<double>(Parameter::MOL_NUM);
  for(int i = 0; i < Parameter::MOL_NUM; i++){
    cm_pos.x += pos[i].x; 
    cm_pos.y += pos[i].y; 
    cm_pos.z += pos[i].z;
  }
  cm_pos.x *= inv_molnum;
  cm_pos.y *= inv_molnum;
  cm_pos.z *= inv_molnum;
  
  double sd_rad = 0.0;
  for(int i = 0; i < Parameter::MOL_NUM; i++){
    const double3 cm2pos = {pos[i].x - cm_pos.x,
			    pos[i].y - cm_pos.y,
			    pos[i].z - cm_pos.z};
    sd_rad += cm2pos.x * cm2pos.x + cm2pos.y * cm2pos.y + cm2pos.z * cm2pos.z;
  }
  sd_rad *= inv_molnum;
  fprintf(fp, "%d %g \n", time, sd_rad);
}

void calculate_tempera(const double3* vel, FILE* fp){
  double sum_vel = 0.0;
  for(int i = 0; i < Parameter::MOL_NUM; i++) 
    sum_vel += vel[i].x * vel[i].x + vel[i].y * vel[i].y + vel[i].z * vel[i].z;
  sum_vel /= 3.0 * Parameter::MOL_NUM;
  fprintf(fp, "%g \n", sum_vel);
}

namespace{
  double3 pos[Parameter::MOL_NUM], vel[Parameter::MOL_NUM], force[Parameter::MOL_NUM];
  Plist pair_list;
}

int main(){
  PRNG prng(1);
  
  const int all_time = 2000000, time_step = 100;
  FILE *fp_rg = fopen("./rg.txt", "w"), *fp_sdrad = fopen("./sdrad.txt", "w"), *fp_temp = fopen("./tempera.txt", "w");
  initialize(pos, vel, force, pair_list);
  create_init_config(pos, vel, prng);
  Mknearlist::create_nearlist(pos, pair_list);
  for(int time = 0; time < all_time; time++){
    update_pos_vel(pos, vel, force);
    if(! Mknearlist::pairlist_is_valid(vel, pair_list ) ){
      std::cerr << "pairlist is re-constructed. \n";
      Mknearlist::create_nearlist(pos, pair_list); 
    }
    clear_force(force);
    F_calculator::calculate_nonbonded_force(pos, force, pair_list);
    F_calculator::calculate_bonded_force(pos, force);
    F_calculator::calculate_random_force(force, prng);
    update_vel(vel, force);

    if(time % time_step == 0){
      calculate_rg(pos, time, fp_rg);
      calculate_dev_rad(pos, time, fp_sdrad);
      calculate_tempera(vel, fp_temp);
    }
  }
}

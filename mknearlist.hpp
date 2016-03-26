#pragma once

#include "parameter.hpp"
#include "mvector3.hpp"

class Mknearlist {
  const double3 eps = {1.0e-5, 1.0e-5, 1.0e-5};
  double3 cell_leng, icell_leng;
  std::array<int, 3> grid_numb;

  void set_grid_info(const double3& L) {
    grid_numb[0] = static_cast<int>(L.x * Parameter::iplist_len);
    grid_numb[1] = static_cast<int>(L.y * Parameter::iplist_len);
    grid_numb[2] = static_cast<int>(L.z * Parameter::iplist_len);
    
    cell_leng.x = L.x / grid_numb[0];
    cell_leng.y = L.y / grid_numb[1];
    cell_leng.z = L.z / grid_numb[2];

    icell_leng = 1.0 / cell_leng;
  }
  
public:
  
  //NOTE: O(N^2) method.
  static void create_nearlist(const double3* pos,
			      Plist& pair_list) {
    for (int i = 0; i < Parameter::MOL_NUM; i++) {
      pair_list[i].clear();
    }
    
    for (int i = 0; i < Parameter::MOL_NUM; i++) {
      const double3 ri = pos[i];
      for (int j = i + 1; j < Parameter::MOL_NUM; j++) {
	const double3 drij = {ri.x - pos[j].x,
			      ri.y - pos[j].y,
			      ri.z - pos[j].z};
	const double dr2 = drij.x * drij.x + drij.y * drij.y + drij.z * drij.z;
	if (dr2 < Parameter::plist_len2) {
	  pair_list[i].push_back(j);
	}
      }
    }
  }
  
  inline int gen_hash(double3& pos, 
		      const double3& origin, 
		      const double3& iL) {
    const double3 origin2pos = pos - origin;
    
    const std::array<int, 3> idx = { static_cast<int>(origin2pos.x - eps.x * icell_leng.x),
				     static_cast<int>(origin2pos.y - eps.y * icell_leng.y),
				     static_cast<int>(origin2pos.z - eps.z * icell_leng.z)};
    return idx[0] + grid_numb[0] * (idx[1] + idx[2] * grid_numb[1]);
  }

  /*  void create_dynamic_celllist(const double3* pos,
			       Plist& pari_list)
  {
    double3 origin(0.0), up_pos(0.0);
    
    origin = pos[0]; up_pos = pos[0];
    for(int i = 0; i < Parameter::MOL_NUM; i++){
      origin.x = (origin.x > pos[i].x) ? pos[i].x : origin.x;
      origin.y = (origin.y > pos[i].y) ? pos[i].y : origin.y;
      origin.z = (origin.z > pos[i].z) ? pos[i].z : origin.z;

      up_pos.x = (up_pos.x < pos[i].x) ? pos[i].x : up_pos.x;
      up_pos.y = (up_pos.y < pos[i].y) ? pos[i].y : up_pos.y;
      up_pos.z = (up_pos.z < pos[i].z) ? pos[i].z : up_pos.z;
    }

    const double3 L = up_pos - origin;
    const double3 iL = 1.0 / L;
    
    set_grid_info(L);
    
    
    }*/

  static bool pairlist_is_valid(const double3* vel,
				const Plist& pair_list) {
    static double max_disp = 0.0;
    double vmax2 = 0.0;
    for (int i = 0; i < Parameter::MOL_NUM; i++) {
      const double temp2 = vel[i].x * vel[i].x + vel[i].y * vel[i].y + vel[i].z * vel[i].z;
      vmax2 = (vmax2 < temp2) ? temp2 : vmax2;
    }
    const double vmax = std::sqrt(vmax2);
  
    max_disp += 2.0 * vmax * Parameter::dt;
    if (max_disp < Parameter::margin) {
      return true;
    } else {
      max_disp = 0.0;
      return false;
    }
  }
};

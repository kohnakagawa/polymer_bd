#pragma once

#include "parameter.hpp"
#include "mvector3.hpp"

class Mknearlist{
public:
  static void create_nearlist(const double3* pos,
		       Plist& pair_list)
  {
    for(int i = 0; i < Parameter::MOL_NUM; i++){
      pair_list[i].clear();
    }
    
    for(int i = 0; i < Parameter::MOL_NUM; i++){
      const double3 ri = pos[i];
      for(int j = i + 1; j < Parameter::MOL_NUM; j++){
	const double3 drij = {ri.x - pos[j].x,
			      ri.y - pos[j].y,
			      ri.z - pos[j].z};
	const double dr2 = drij.x * drij.x + drij.y * drij.y + drij.z * drij.z;
	if(dr2 < Parameter::plist_len2){
	  pair_list[i].push_back(j);
	}
      }
    }
  }

  static bool pairlist_is_valid(const double3* vel,
			 const Plist& pair_list)
  {
    static double max_disp = 0.0;
    double vmax2 = 0.0;
    for(int i = 0; i < Parameter::MOL_NUM; i++){
      const double temp2 = vel[i].x * vel[i].x + vel[i].y * vel[i].y + vel[i].z * vel[i].z;
      vmax2 = (vmax2 < temp2) ? temp2 : vmax2;
    }
    const double vmax = std::sqrt(vmax2);
  
    max_disp += 2.0 * vmax * Parameter::dt;
    if(max_disp < Parameter::margin){
      return true;
    }else{
      max_disp = 0.0;
      return false;
    }
  }
  
};

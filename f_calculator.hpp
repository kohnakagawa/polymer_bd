#pragma once

#include "parameter.hpp"
#include "mvector3.hpp"
#include "prng.hpp"

class F_calculator{
public:
  static void calculate_bonded_force(const double3* __restrict pos,
				     double3* __restrict force)
  {
    for(int i = 0; i < Parameter::MOL_NUM - 1; i++){
      const double3 dr		= pos[i + 1] - pos[i];
      const double  dr2		= dr.x * dr.x + dr.y * dr.y + dr.z * dr.z;
      const double  inv_dr	= 1.0 / std::sqrt(dr2);
      const double  cf_bond	= Parameter::cf_b * (inv_dr - Parameter::inv_b);
      const double3 Fbond(cf_bond * dr.x, cf_bond * dr.y, cf_bond * dr.z);
    
      force[i + 1] += Fbond; 
      force[i    ] -= Fbond;
    }
  }
  
  static void calculate_nonbonded_force(const double3* __restrict pos,
					double3* __restrict force,
					const Plist& pair_list)
  {
    for(int i = 0; i < Parameter::MOL_NUM; i++){
      const int pair_num = pair_list[i].size();
      const double3 ri = pos[i];
      for(int j = 0; j < pair_num; j++){
	const int jid = pair_list[i][j];
	const double3 drij = {ri.x - pos[jid].x, 
			      ri.y - pos[jid].y,
			      ri.z - pos[jid].z};
	const double dr2 = drij.x * drij.x + drij.y * drij.y + drij.z * drij.z;
	if(dr2 < Parameter::cutof2){
	  const double inv_dr2 = 1.0 / dr2;
	  const double inv_dr6 = inv_dr2 * inv_dr2 * inv_dr2;
	
	  const double dF_norm = Parameter::cf_nbnded * inv_dr6 * 6.0 * (2.0 * inv_dr6 - 1.0) * inv_dr2;
	  const double3 dF_nb(dF_norm * drij.x, dF_norm * drij.y, dF_norm * drij.z);
	  force[i] += dF_nb;
	  force[jid] -= dF_nb;
	}
      }
    }
  }

  static void calculate_random_force(double3* force,
				     PRNG& prng)
  {
    for(int i = 0; i < Parameter::MOL_NUM; i++){
      const double3 Rforce(prng.Normal(0) * Parameter::nois_amp,
			   prng.Normal(0) * Parameter::nois_amp,
			   prng.Normal(0) * Parameter::nois_amp);
      force[i] += Rforce;
    }
  }
};

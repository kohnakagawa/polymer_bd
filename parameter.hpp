#pragma once

#include <array>
#include <vector>
#include <cmath>

class Parameter {
public:
  enum{MOL_NUM = 100,};  

  static constexpr double bond = 1.0, inv_b  = 1.0 / bond;
  static constexpr double rc   = 1.0, inv_rc = 1.0 / rc;
  static constexpr double cf_b = 50.0;
  static constexpr double dt   = 0.0005;
  static constexpr double gamm = 0.05;
  static constexpr double mass = 1.0;

  static constexpr double time_cf  = 1.0 / (1.0 + 0.5 * gamm * dt);
  static constexpr double tempera  = 1.0;
  static constexpr double half_dt_mass  = 0.5 * dt / mass;
  static constexpr double nois_amp = std::sqrt(2.0 * tempera * gamm * mass / dt);
  
  static constexpr double cutof = 1.5;
  static constexpr double cutof2 = cutof * cutof;
  static constexpr double cf_nbnded = 1.0;

  static constexpr double plist_len  = 1.8;
  static constexpr double plist_len2 = plist_len * plist_len;
  static constexpr double iplist_len = 1.0 / plist_len;
  static constexpr double margin = plist_len - cutof;
};

typedef std::array<std::vector<int>, Parameter::MOL_NUM> Plist;

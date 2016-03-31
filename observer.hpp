#pragma once

#include <string>
#include <cstdlib>
#include "mvector3.hpp"

class Observer {
  enum {
    RadGyr = 0,
    SdRad,
    Tempera,
    Config,
    Num_Files,
  };

  FILE* fp[Num_Files] = {nullptr};
  std::string fid_to_name(const int fid) const ;
  
public:
  Observer() {
    for (int i = 0; i < Num_Files; i++) {
      fp[i] = fopen(fid_to_name(i).c_str(), "w");
    }
  }
  ~Observer() {
    for (int i = 0; i < Num_Files; i++) {
      fclose(fp[i]);
    }
  }
  
  void calculate_rg(const double3* pos, const int time);
  void calculate_dev_rad(const double3* pos, const int time);
  void calculate_tempera(const double3* vel);
  void print_config(const double3* pos, const int time);
};

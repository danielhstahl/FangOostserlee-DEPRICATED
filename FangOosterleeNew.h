#ifndef __FANGOOSTERLEENEW_H_INCLUDED__
#define __FANGOOSTERLEENEW_H_INCLUDED__
#define _USE_MATH_DEFINES
#include <cmath>
#include "Complex.h"
#include <vector>

class FangOost{
  double xmin;
  double xmax;
  int h;//discrete x
  int k;//discrete u
  double du;
  double dx;
  std::vector<double> f;
  public:
    FangOost(double, double, int, int);
    void set_values(double, double, int, int);
    void computeInv(
      auto&& fnInv //function to invert...function of u only
    );
    std::vector<double> computeConvolution(
      auto& //vk as defined in fangoosterlee
    );
    void computeExpectation(
      auto& //callback takes (f, du, dx, k, h)
    );
};

#include "FangOosterleeNew.hpp"
#endif

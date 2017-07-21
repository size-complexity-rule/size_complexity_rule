#ifndef __RATES_HPP__
#define __RATES_HPP__

#include <cmath>

#include <vector>
using std::vector;

#include "population.hpp"

#include "random.hpp"

enum ReactionType { aggregation=0, reproduction, death, dissociation };

class Rates {
protected:
  vector<vector<double>> aggregation;
  vector<double> reproduction;
  vector<double> death;
  vector<double> dissociation;

  std::vector<double> fitness;

  double total_rate;

  double k_reproduction, k_aggregation, k_death, k_dissociation;
  double stickiness, C;

  double R1(double C, double p) {
    double x = 1.;
    for (int i = 0; i < 10; i++)
      x = x - (exp((1.-p)*x*x)*p*p + 4.*C*x*(-3. + 2.*(1.-p)*x*x))/(-2.*((p-1.)*x*exp((1.-p)*x*x)*p*p + 6.*C*(1. + 2.*(p-1.)*x*x)));
    return x*x*(1. - 8.*C/(p*p)*x*exp((p-1.)*x*x));
  }

  int lmax;
  int maxl2;
public:
  Rates (double _k_reproduction, double _k_aggregation, double _k_death,
         double _k_dissociation, double _stickiness, double _C, int p) {
    k_reproduction = _k_reproduction;
    k_aggregation  = _k_aggregation;
    k_death        = _k_death;
    k_dissociation = _k_dissociation;

    stickiness     = _stickiness;

    lmax = 1;

    // Allocate the vectors
    reproduction  .resize(max_size, 0.0);
    death         .resize(max_size, 0.0);
    dissociation  .resize(max_size, 0.0);
    aggregation   .resize(max_size);
    for (auto& element: aggregation) element.resize(max_size, 0.0);

    // Define the fitnesses
    C = _C;
    fitness.resize(max_size);
    fitness[0] = 0.0;
    fitness[1] = R1(C, p);
    for (int i = 2; i < p; i++)
      fitness[i] = double(i) * fitness[1];
    for (int i = p; i < max_size; i++) {
      double k = i%p;
      double m = i/p;
      fitness[i] = double(i) * pow(double(p),4.)*pow(m+1.,4.*k/p)*pow(m,4.*(p-k)/p)/(pow(double(i),4.)*432.*C*C);
    }
  }

  virtual void update(Population* pop) =0;

  double time_to_next_reaction () {
    if (total_rate > 0.0) {
      return -log(uniform())/total_rate;
    } else {
      exit(0);
    }
  }

  void choose_reaction(ReactionType& reaction, int& l1, int& l2) {
    double threshold = uniform()*total_rate;
    double rate = 0.0;

    reaction = ReactionType::aggregation;
    l2 = 1;
    reaction = ReactionType::aggregation;
    for (l1 = 1; l1 <=lmax; l1++)
      for (l2 = 1; l2 <= l1 and l2 <= maxl2; l2++) {
        rate += aggregation[l1][l2];
        if (rate > threshold) return;
      }
    l2 = 0;


    reaction = ReactionType::reproduction;
    for (l1 = 1; l1 <= lmax; l1++) {
      rate += reproduction[l1];
      if (rate > threshold) return;
    }

    reaction = ReactionType::death;
    for (l1 = 1; l1 <= lmax; l1++) {
      rate += death[l1];
      if (rate > threshold) return;
    }

    reaction = ReactionType::dissociation;
    for (l1 = 2; l1 <= lmax; l1++) {
      rate += dissociation[l1];
      if (rate > threshold) return;
    }
  }
};


class SphericalRates: public Rates {
public:
  SphericalRates(double _k_reproduction, double _k_aggregation, double _k_death,
                 double _k_dissociation, double _stickiness, double _C, int p)
  : Rates(_k_reproduction, _k_aggregation, _k_death, _k_dissociation, _stickiness, _C, p) {
    maxl2 = 1;
  }

  void update(Population* pop) {
    total_rate = 0.0;
    lmax = pop->get_lmax();

    // Aggregation
    int l2 = 1;
    for (int l1 = 1; l1 <= lmax; l1++) {
      double rate;
      if (l1 == 1) {
        rate = k_aggregation*pop->size(1)*(pop->size(1)-1)*stickiness*stickiness*pow(double(l1), 2.0/3.0);
      } else {
        rate = k_aggregation*pop->size(l1)*pop->size(1)*stickiness*stickiness*pow(double(l1), 2.0/3.0);
      }
      aggregation[l1][l2] = rate;
      total_rate += rate;
    }

    // Reproduction
    for (size_t l = 1; l <= lmax; l++) {
      double rate = k_reproduction*fitness[l]*pop->size(l);
      reproduction[l] = rate;
      total_rate += rate;
    }

    // Death
    for (size_t l = 1; l <= lmax; l++) {
      double rate = k_death*l*pop->size(l)*pop->total_size();
      death[l] = rate;
      total_rate += rate;
    }

    // Dissociation
    for (size_t l = 2; l <= lmax; l++) {
      double rate = k_dissociation*pop->size(l)*pow(l,(2.0/3.0));
      dissociation[l] = rate;
      total_rate += rate;
    }
  }
};

class LinearRates: public Rates {
public:
  LinearRates(double _k_reproduction, double _k_aggregation, double _k_death,
              double _k_dissociation, double _stickiness,    double _C, int p)
  : Rates(_k_reproduction, _k_aggregation, _k_death, _k_dissociation, _stickiness, _C, p) {
    maxl2 = max_size;
  }

  void update(Population* pop) {
    total_rate = 0.0;
    lmax = pop->get_lmax();

    // Aggregation
    for (size_t l1 = 1; l1 <= lmax; l1++)
      for (size_t l2 = 1; l2 <= l1; l2++) {
        double rate;
        if (l1 != l2) {
          rate = k_aggregation*pop->size(l1)*pop->size(l2)*stickiness*stickiness;
        } else {
          rate = k_aggregation*pop->size(l1)*(pop->size(l1)-1)*stickiness*stickiness;
        }
        aggregation[l1][l2] = rate;
        total_rate += rate;
      }

    // Reproduction
    for (size_t l = 1; l <= lmax; l++) {
      double rate = k_reproduction*fitness[l]*pop->size(l);
      reproduction[l] = rate;
      total_rate += rate;
    }

    // Death
    for (size_t l = 1; l <= lmax; l++) {
      double rate = k_death*l*pop->size(l)*pop->total_size();
      death[l] = rate;
      total_rate += rate;
    }

    // Dissociation
    for (size_t l = 2; l <= lmax; l++) {
      double rate = k_dissociation*(l-1)*pop->size(l);
      dissociation[l] = rate;
      total_rate += rate;
    }
  }
};

#endif /* end of include guard: __RATES_HPP__ */

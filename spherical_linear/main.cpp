
#include <iostream>
#include <fstream>
#include <functional>

const size_t max_size = 65536;

enum Geometry { Spherical=0, Linear=1 };

#include "random.hpp"
#include "utilities.hpp"

#include "rates.hpp"
#include "population.hpp"

using namespace std;

void save_results(Population* pop, int iterate, double t, std::ofstream& data, ReactionType reaction, int l);
bool save(int iterate);

void update_spherical(Population* pop, ReactionType next_reaction, int l1, int l2);
void update_linear(Population* pop, ReactionType next_reaction, int l1, int l2);

int main(int argc, char const *argv[]) {
  initialize_random();

  //////////////////////////////////////////////////////////////////////////////
  //   Read parameters
  if (argc != 11) {
    std::cout << "error: wrong options!" << std::endl;
    exit(1);
  }

  double k_agg              = atof(argv[2]);
  double k_dis              = atof(argv[3]);
  double k_rep              = atof(argv[4]);
  double k_dea              = atof(argv[5]);
  double stickiness         = atof(argv[6]);
  double C                  = atof(argv[7]);
  int    p                  = atoi(argv[10]);

  Geometry geometry;
  string geometry_name = string(argv[1]);
  if (geometry_name == "linear") {
    geometry = Geometry::Linear;
  } else if (geometry_name == "spherical") {
    geometry = Geometry::Spherical;
  } else {
    std::cout << "error: unrecognized geometry type (" << geometry_name << ")" << std::endl;
    exit(0);
  }

  size_t initial_population = atoi(argv[8]);

  size_t max_iterate        = atoi(argv[9]);

  double probability_single = 1.0 - stickiness*stickiness;
  //////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////
  //   Files
  string geometry_type;
  switch (geometry) {
    case Geometry::Linear:
      geometry_type = "linear";
      break;
    case Geometry::Spherical:
      geometry_type = "spherical";
      break;
  }

  string filename
          = string("_p")  + to_string(p)
          + string("_kp") + to_string(k_agg)
          + string("_km") + to_string(k_dis)
          + string("_kr") + to_string(k_rep)
          + string("_kd") + to_string(k_dea)
          + string("_st") + to_string(stickiness)
          + string(".dat");

  std::ofstream data_file(string("results/data_") + geometry_type + filename);

  string stamp
          = string("#"                  )                                + string("\n")
          + string("# Run at "          ) + get_time_and_date()          + string("\n")
          + string("#     on "          ) + get_hostname()               + string("\n")
          + string("#"                  )                                + string("\n")
          + string("# Parameters:"      )                                + string("\n")
          + string("#     k+          " ) + to_string(k_agg)             + string("\n")
          + string("#     k-          " ) + to_string(k_dis)             + string("\n")
          + string("#     kR          " ) + to_string(k_rep)             + string("\n")
          + string("#     kD          " ) + to_string(k_dea)             + string("\n")
          + string("#     stickiness  " ) + to_string(stickiness)        + string("\n")
          + string("#     p_single    " ) + to_string(probability_single)+ string("\n")
          + string("#     p           " ) + to_string(p)                 + string("\n");

  data_file
          << "##################################################################\n"
          << "# General data file for " << geometry_type << " group\n"
          << stamp
          << "##################################################################\n\n"
          << "# iterate\ttime\tlmax\taverage_size\tgroup_number\ttotal_cells_number\treaction\tl"
          << std::endl;
  //////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////
  //   Initialization
  Rates* rates;
  function<void(Population*, ReactionType, size_t, size_t)> update;
  switch (geometry) {
    case Geometry::Linear:
      rates  = new LinearRates(k_rep, k_agg, k_dea, k_dis, stickiness, C, p);
      update = update_linear;
      break;
    case Geometry::Spherical:
      rates  = new SphericalRates(k_rep, k_agg, k_dea, k_dis, stickiness, C, p);
      update = update_spherical;
      break;
  }

  Population* population = new Population(probability_single);
  for (int i = 0; i < initial_population; i++) population->add(1);
  //////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////
  //   Main cycle
  double t = 0.0;
  for (size_t t1 = 1; t1 <= max_iterate; t1++) {
    rates->update(population);

    t += rates->time_to_next_reaction();

    int l1, l2;
    ReactionType next_reaction;
    rates->choose_reaction(next_reaction, l1, l2);
    update(population, next_reaction, l1, l2);

    if ( save(t1) ) save_results(population, t1, t, data_file, next_reaction, l1);
  }
  //////////////////////////////////////////////////////////////////////////////
  
  data_file << flush;

  return 0;
}

void save_results(Population* pop, int iterate, double t, std::ofstream& data, ReactionType reaction, int l) {
  const  string sep = string("\t");
  double average_size = double(pop->total_size())/pop->size();
  int lmax = pop->get_lmax();

  data << iterate << sep << t << sep
       << lmax                << sep
       << average_size        << sep
       << pop->size()         << sep
       << pop->total_size()   << sep
       << reaction            << sep
       << l                   << "\n";
}

bool save(int iterate) {
  return !(iterate%1000);
}

void update_spherical(Population* pop, ReactionType next_reaction, int l1, int l2) {
  switch (next_reaction) {
    case ReactionType::aggregation: {
      pop->remove(l1);
      pop->remove(l2);
      pop->add(l1+l2);
    } break;
    case ReactionType::reproduction: {
      if (uniform() < pop->get_split_probability()) {
        pop->add(1);
      } else {
        pop->remove(l1);
        pop->add(l1+1);
      }
    } break;
    case ReactionType::death: {
      pop->remove(l1);
      if (l1 > 1) pop->add(l1-1);
    } break;
    case ReactionType::dissociation: {
      pop->remove(l1);
      pop->add(l1-1);
      pop->add(1);
    } break;
  }

  pop->update_lmax();
}

void update_linear(Population* pop, ReactionType next_reaction, int l1, int l2) {
  switch (next_reaction) {
    case ReactionType::aggregation: {
      pop->remove(l1);
      pop->remove(l2);
      pop->add(l1+l2);
    } break;
    case ReactionType::reproduction: {
      if (uniform() < pop->get_split_probability()) {
        pop->add(1);
      } else {
        pop->remove(l1);
        pop->add(l1+1);
      }
    } break;
    case ReactionType::death: {
      pop->remove(l1);

      int dead = random_integer(1, l1);
      if (l1 - dead > 0) {
        pop->add(l1-dead);
      }
      if (dead - 1 > 0) {
        pop->add(dead-1);
      }
    } break;
    case ReactionType::dissociation: {
      if (l1 > 1) {
        size_t break_point = random_integer(1, l1-1);
        pop->remove(l1);
        pop->add(break_point);
        pop->add(l1-break_point);
      }
    } break;
  }

  pop->update_lmax();
}

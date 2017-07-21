#ifndef __POPULATION_HPP__
#define __POPULATION_HPP__

#include <vector>
#include <iostream>

class Population {
private:
  std::vector<int> pop;
  int pop_size;
  int lmax;
  int total_cells;

  double split_probability;
public:
  Population (double _split_probability) {
    split_probability = _split_probability;

    lmax = 1;
    pop_size = 0;
    total_cells = 0;
    pop.resize(max_size, 0);
  }

  void add(int l) {
    pop[l]   ++;
    pop_size ++;
    total_cells += l;
  }

  void remove(int l) {
    pop[l]   --;
    pop_size --;
    total_cells -= l;
  }

  int size() {
    return pop_size;
  }

  int total_size() {
    return total_cells;
  }

  int size(size_t i) {
    return pop[i];
  }

  int get_lmax() {
    return lmax;
  }

  double get_split_probability() {
    return split_probability;
  }

  void print(double t=0) {
    std::cout << t << "\t" << pop_size << "\t" << lmax << "\t" << double(total_cells)/pop_size << "\t" << total_cells << std::endl;
  }

  void update_lmax() {
    int limit = 2*(lmax+1);
    if(max_size < limit) limit = max_size;
    lmax = 0;
    for (int l = 1; l < limit; l++)
      if (pop[l] > 0) lmax = l;
    lmax ++;
  }
};

#endif /* end of include guard: __POPULATION_HPP__ */

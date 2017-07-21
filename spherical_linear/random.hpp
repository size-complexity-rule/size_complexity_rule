#ifndef __RANDOM_H__
#define __RANDOM_H__

#include <random>

void initialize_random(int test = 0);

double  uniform();
int     random_integer(int min, int max);

using generator_t = std::ranlux48;
generator_t generator;
std::uniform_real_distribution<double>  uniform_distribution;

double uniform() {
    return uniform_distribution(generator);
}

int random_integer(int min, int max) {
    std::uniform_int_distribution<int> distribution(min, max);
    return distribution(generator);
}

int random_integer(int max) {
    std::uniform_int_distribution<int> distribution(0, max);
    return distribution(generator);
}

void initialize_random(int test) {
    std::random_device rd;
    if (not test) {
      generator = generator_t( rd() );
    } else {
      generator = generator_t( 213158486757837 );
    }

    uniform_distribution = std::uniform_real_distribution<double>(0, 1);
}

#endif /* end of include guard: __RANDOM_H__ */

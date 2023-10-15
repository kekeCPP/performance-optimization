/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "vector.hpp"
#include <vector>

#if !defined(ANALYSIS_HPP)
#define ANALYSIS_HPP

namespace Analysis {
std::vector<double> correlation_coefficients_p(std::vector<Vector> datasets, unsigned int MAX_THREADS);
std::vector<double> correlation_coefficients(std::vector<Vector> datasets);
double pearson(Vector vec1, Vector vec2);
};

#endif

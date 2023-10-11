/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "filters.hpp"
#include "matrix.hpp"
#include "ppm.hpp"
#include <cstdlib>
#include <iostream>

int main(int argc, char const* argv[])
{
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " [infile] [outfile] [num_threads]" << std::endl;
        std::exit(1);
    }

    PPM::Reader reader {};
    PPM::Writer writer {};

    auto m { reader(argv[1]) };
    Filter::threshold(m, std::stoul(15));

    writer(m, argv[2]);

    return 0;
}
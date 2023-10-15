/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "analysis.hpp"
#include "dataset.hpp"
#include <iostream>
#include <cstdlib>
#include <chrono>

int main(int argc, char const* argv[])
{
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " [dataset] [outfile] [number_of_threads]" << std::endl;
        std::exit(1);
    }

    auto datasets { Dataset::read(argv[1]) };

    auto corrs { Analysis::correlation_coefficients_p(datasets, strtol(argv[3], NULL, 10)) };

    // auto corrs2 { Analysis::correlation_coefficients(datasets) };

    // for (int i = 0; i < corrs.size(); i++) {
    //     std::cout << "PAR[" << i << "]: " << corrs[i] << "\nNON_PAR[" << i << "]: " << corrs2[i] << "\n";
    //     if (corrs[i] != corrs2[i]) { std::cout << "FALSE!!!\n"; }
    // }
    
    auto start = std::chrono::high_resolution_clock::now();
    Dataset::write_p(corrs, argv[2]);
    auto stop = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout << "TIME TAKEN: " << duration.count() << "\n";

    return 0;
}

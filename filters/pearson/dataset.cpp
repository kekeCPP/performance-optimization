/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "dataset.hpp"
#include "vector.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <iomanip>
#include <limits>
#include <pthread.h>

namespace Dataset {
std::vector<Vector> read(std::string filename)
{
    unsigned dimension {};
    std::vector<Vector> result {};
    std::ifstream f {};

    f.open(filename);

    if (!f) {
        std::cerr << "Failed to read dataset(s) from file " << filename << std::endl;
        return result;
    }

    f >> dimension;
    std::string line {};

    std::getline(f, line); // ignore first newline

    while (std::getline(f, line)) {
        std::stringstream ss { line };
        Vector new_vec { dimension };
        std::copy_n(std::istream_iterator<double> { ss },
            dimension,
            new_vec.get_data());
        result.push_back(new_vec);
    }

    return result;
}

void write(std::vector<double> data, std::string filename)
{
    std::ofstream f {};

    f.open(filename);

    if (!f) {
        std::cerr << "Failed to write data to file " << filename << std::endl;
        return;
    }

    for (auto i { 0 }; i < data.size(); i++) {
        f << std::setprecision(std::numeric_limits<double>::digits10 + 1) << data[i] << std::endl;
    }
}

struct thread_data {
    unsigned int thread_id;
    std::vector<double>* data;
    std::string filename;
    std::ofstream* f;
    unsigned int start_index;
    unsigned int size;
};

void* write_par(void* thread_args)
{
    thread_data* my_data;
    my_data = (thread_data*) thread_args;

    (*my_data->f).open(my_data->filename);

    if (!(*my_data->f)) {
        std::cerr << "Failed to write data to file " << my_data->filename << std::endl;
    }
    else {
        // unsigned int index = 0;
        for (int i = my_data->start_index; i < (my_data->start_index + my_data->size); i++) {
            (*my_data->f).seekp(i);
            (*my_data->f) << std::setprecision(std::numeric_limits<double>::digits10 + 1) << (*my_data->data)[i] << "\n";
        }
    }

}

void write_p(std::vector<double> data, std::string filename)
{
    int MAX_THREADS = 4;
    thread_data* thread_data_array = new thread_data[MAX_THREADS];
    pthread_t* p_threads = new pthread_t[MAX_THREADS];
    std::ofstream f {};

    for (int i = 0; i < MAX_THREADS; i++) {
        thread_data_array[i].thread_id = i;
        thread_data_array[i].data = &data;
        thread_data_array[i].filename = filename;
        thread_data_array[i].f = &f;
        thread_data_array[i].start_index = (data.size() / MAX_THREADS) * i;
        thread_data_array[i].size = data.size() / MAX_THREADS;

        pthread_create(&p_threads[i], NULL, write_par, (void*) &thread_data_array[i]);
    }

    for (int i = 0; i < MAX_THREADS; i++) {
        pthread_join(p_threads[i], NULL);
    }

    delete[] thread_data_array;
    delete[] p_threads;
}

};

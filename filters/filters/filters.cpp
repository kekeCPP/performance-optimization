/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "filters.hpp"
#include "matrix.hpp"
#include "ppm.hpp"
#include <cmath>

namespace Filter {

namespace Gauss {
    void get_weights(int n, double* weights_out)
    {
        for (auto i { 0 }; i <= n; i++) {
            double x { static_cast<double>(i) * max_x / n };
            weights_out[i] = exp(-x * x * pi);
        }
    }
}

Matrix blur(Matrix &dst, const int radius)
{
    Matrix scratch { PPM::max_dimension };

    double w[Gauss::max_radius] {};
    Gauss::get_weights(radius, w);

    //cache value frequently used, never changed
    const auto dstXsize = dst.get_x_size();
    const auto dstYSize = dst.get_y_size();

    //pointers for r,g,b in dst matrix
    auto dstR = dst.get_R();
    auto dstG = dst.get_G();
    auto dstB = dst.get_B();

    //non constant pointers so values can be changed
    auto dstRnonCon = dst.get_R_nonconst();
    auto dstGnonCon = dst.get_G_nonconst();
    auto dstBnonCon = dst.get_B_nonconst();

    //pointers for r,g,b scratch matrix
    auto scrR = scratch.get_R();
    auto scrG = scratch.get_G();
    auto scrB = scratch.get_B();

    //non constant pointers so values can be changed
    auto scrRnonCon = scratch.get_R_nonconst();
    auto scrGnonCon = scratch.get_G_nonconst();
    auto scrBnonCon = scratch.get_B_nonconst();

    const auto scrXsize = scratch.get_x_size();

    for (auto x { 0 }; x < dstXsize; x++) {
        for (auto y { 0 }; y < dstYSize; y++) {
            auto r { w[0] * dst.r(x, y) }, g { w[0] * dst.g(x, y) }, b { w[0] * dst.b(x, y) }, n { w[0] };

            for (auto wi { 1 }; wi <= radius; wi++) {
                auto wc { w[wi] };
                auto x2 { x - wi };
                if (x2 >= 0) {
                    r += wc * dstR[y * dstXsize + x2];
                    g += wc * dstG[y * dstXsize + x2];
                    b += wc * dstB[y * dstXsize + x2];
                    n += wc;
                }
                x2 = x + wi;
                if (x2 < dstXsize) {
                    r += wc * dstR[y * dstXsize + x2];
                    g += wc * dstG[y * dstXsize + x2];
                    b += wc * dstB[y * dstXsize + x2];
                    n += wc;
                }
            }
            scrRnonCon[y * scrXsize + x] = r / n;
            scrGnonCon[y * scrXsize + x] = g / n;
            scrBnonCon[y * scrXsize + x] = b / n;
        }
    }

    for (auto x { 0 }; x < dstXsize; x++) {
        for (auto y { 0 }; y < dstYSize; y++) {
            auto r { w[0] * scratch.r(x, y) }, g { w[0] * scratch.g(x, y) }, b { w[0] * scratch.b(x, y) }, n { w[0] };

            for (auto wi { 1 }; wi <= radius; wi++) {
                auto wc { w[wi] };
                auto y2 { y - wi };
                if (y2 >= 0) {
                    r += wc * scrR[y2 * scrXsize + x];
                    g += wc * scrG[y2 * scrXsize + x];
                    b += wc * scrB[y2 * scrXsize + x];
                    n += wc;
                }
                y2 = y + wi;
                if (y2 < dstYSize) {
                    r += wc * scrR[y2 * scrXsize + x];
                    g += wc * scrG[y2 * scrXsize + x];
                    b += wc * scrB[y2 * scrXsize + x];
                    n += wc;
                }
            }

            dstRnonCon[y * dstXsize + x] = r / n;
            dstGnonCon[y * dstXsize + x] = g / n;
            dstBnonCon[y * dstXsize + x] = b / n;
        }
    }

    return dst;
}

Matrix threshold(Matrix &m)
{
    unsigned sum {}, nump { m.get_x_size() * m.get_y_size() };

    //pointers for r,g,b in dst matrix
    auto dstR = m.get_R();
    auto dstG = m.get_G();
    auto dstB = m.get_B();

    for (auto i { 0 }; i < nump; i+=4) {
        sum += dstR[i] + dstG[i] + dstB[i];
        sum += dstR[i + 1] + dstG[i + 1] + dstB[i + 1];
        sum += dstR[i + 2] + dstG[i + 2] + dstB[i + 2];
        sum += dstR[i + 3] + dstG[i + 3] + dstB[i + 3];
    }

    sum /= nump;

    unsigned psum {};

    //non constant pointers so values can be changed
    auto dstR2 = m.get_R_nonconst();
    auto dstG2 = m.get_G_nonconst();
    auto dstB2 = m.get_B_nonconst();

    for (auto i { 0 }; i < nump; i++) {
        //psum = dst.r(i, 0) + dst.g(i, 0) + dst.b(i, 0);
        psum = dstR[i] + dstG[i] + dstB[i];
        if (sum > psum) {
            //dst.r(i, 0) = dst.g(i, 0) = dst.b(i, 0) = 0;
            dstR2[i] = dstG2[i] = dstB2[i] = 0;
        } else {
            //dst.r(i, 0) = dst.g(i, 0) = dst.b(i, 0) = 255;
            dstR2[i] = dstG2[i] = dstB2[i] = 255;
        }
    }

    return 0;
}

//parallelised versions
Matrix blur_par(Matrix &dst, const int radius, const int MAX_THREADS)
{
    std::cout << "created: " << MAX_THREADS << " threads\n";

    Matrix scratch { PPM::max_dimension };

    double w[Gauss::max_radius] {};
    Gauss::get_weights(radius, w);

    //cache value frequently used, never changed
    const auto dstXsize = dst.get_x_size();
    const auto dstYSize = dst.get_y_size();

    //pointers for r,g,b in dst matrix
    auto dstR = dst.get_R();
    auto dstG = dst.get_G();
    auto dstB = dst.get_B();

    //non constant pointers so values can be changed
    auto dstRnonCon = dst.get_R_nonconst();
    auto dstGnonCon = dst.get_G_nonconst();
    auto dstBnonCon = dst.get_B_nonconst();

    //pointers for r,g,b scratch matrix
    auto scrR = scratch.get_R();
    auto scrG = scratch.get_G();
    auto scrB = scratch.get_B();

    //non constant pointers so values can be changed
    auto scrRnonCon = scratch.get_R_nonconst();
    auto scrGnonCon = scratch.get_G_nonconst();
    auto scrBnonCon = scratch.get_B_nonconst();

    const auto scrXsize = scratch.get_x_size();

    for (auto x { 0 }; x < dstXsize; x++) {
        for (auto y { 0 }; y < dstYSize; y++) {
            auto r { w[0] * dst.r(x, y) }, g { w[0] * dst.g(x, y) }, b { w[0] * dst.b(x, y) }, n { w[0] };

            for (auto wi { 1 }; wi <= radius; wi++) {
                auto wc { w[wi] };
                auto x2 { x - wi };
                if (x2 >= 0) {
                    r += wc * dstR[y * dstXsize + x2];
                    g += wc * dstG[y * dstXsize + x2];
                    b += wc * dstB[y * dstXsize + x2];
                    n += wc;
                }
                x2 = x + wi;
                if (x2 < dstXsize) {
                    r += wc * dstR[y * dstXsize + x2];
                    g += wc * dstG[y * dstXsize + x2];
                    b += wc * dstB[y * dstXsize + x2];
                    n += wc;
                }
            }
            scrRnonCon[y * scrXsize + x] = r / n;
            scrGnonCon[y * scrXsize + x] = g / n;
            scrBnonCon[y * scrXsize + x] = b / n;
        }
    }

    for (auto x { 0 }; x < dstXsize; x++) {
        for (auto y { 0 }; y < dstYSize; y++) {
            auto r { w[0] * scratch.r(x, y) }, g { w[0] * scratch.g(x, y) }, b { w[0] * scratch.b(x, y) }, n { w[0] };

            for (auto wi { 1 }; wi <= radius; wi++) {
                auto wc { w[wi] };
                auto y2 { y - wi };
                if (y2 >= 0) {
                    r += wc * scrR[y2 * scrXsize + x];
                    g += wc * scrG[y2 * scrXsize + x];
                    b += wc * scrB[y2 * scrXsize + x];
                    n += wc;
                }
                y2 = y + wi;
                if (y2 < dstYSize) {
                    r += wc * scrR[y2 * scrXsize + x];
                    g += wc * scrG[y2 * scrXsize + x];
                    b += wc * scrB[y2 * scrXsize + x];
                    n += wc;
                }
            }

            dstRnonCon[y * dstXsize + x] = r / n;
            dstGnonCon[y * dstXsize + x] = g / n;
            dstBnonCon[y * dstXsize + x] = b / n;
        }
    }

    return dst;
}

struct thread_data{
        int thread_id;
        int thread_number;
        int nump;
        const unsigned char* dstR;
        const unsigned char* dstG;
        const unsigned char* dstB;
        int* sum;
};

void *threadFunc(void * thread_arg){
    struct thread_data *my_data;
    my_data = (struct thread_data *) thread_arg;
    int sum = 0;
    for (auto i { my_data->thread_id }; i < my_data->nump; i += my_data->thread_number) {
        *my_data->sum +=  my_data->dstR[i] + my_data->dstG[i] + my_data->dstG[i];
        sum += my_data->dstR[i] + my_data->dstG[i] + my_data->dstG[i];
    }
    //std::cout << *my_data->sum << " ";
    std::cout << sum << " ";

    pthread_exit(NULL);
}

Matrix threshold_par(Matrix &m, const int MAX_THREADS)
{
    unsigned sum {}, nump { m.get_x_size() * m.get_y_size() };

    //pointers for r,g,b in dst matrix
    auto dstR = m.get_R();
    auto dstG = m.get_G();
    auto dstB = m.get_B();

    struct thread_data thread_data_array[MAX_THREADS];
    pthread_t p_threads[MAX_THREADS];
    int thread_sum[MAX_THREADS];

    for(auto i { 0 }; i < MAX_THREADS; i++){
        thread_data_array[i].thread_id = i;
        thread_data_array[i].thread_number = MAX_THREADS;
        thread_data_array[i].nump = nump;
        thread_data_array[i].dstR = dstR;
        thread_data_array[i].dstG = dstG;
        thread_data_array[i].sum = &thread_sum[i];

        pthread_create(
            &p_threads[i],
            NULL,
            threadFunc,
            (void*) &thread_data_array[i]
        );
    }

    int sumReal = 0;
    for (auto i { 0 }; i < nump; i++) {
        sumReal += dstR[i] + dstG[i] + dstB[i];
    }

    for (auto i { 0 } ; i < MAX_THREADS; i++) {
        pthread_join(p_threads[i], NULL); // Wait for all threads to terminate
    }
    std::cout << "\n\n";
    for(auto i = 0; i < MAX_THREADS; i++){
        std::cout << thread_sum[i] << " ";
        sum+= thread_sum[i];
    }
    std::cout << "\n\nsum: "<< sum << "\n\n";

    std::cout << "\n\nReal sum: "<< sumReal << "\n\n";
    sum /= nump;

    unsigned psum {};

    //non constant pointers so values can be changed
    auto dstR2 = m.get_R_nonconst();
    auto dstG2 = m.get_G_nonconst();
    auto dstB2 = m.get_B_nonconst();

    for (auto i { 0 }; i < nump; i++) {
        //psum = dst.r(i, 0) + dst.g(i, 0) + dst.b(i, 0);
        psum = dstR[i] + dstG[i] + dstB[i];
        if (sum > psum) {
            //dst.r(i, 0) = dst.g(i, 0) = dst.b(i, 0) = 0;
            dstR2[i] = dstG2[i] = dstB2[i] = 0;
        } else {
            //dst.r(i, 0) = dst.g(i, 0) = dst.b(i, 0) = 255;
            dstR2[i] = dstG2[i] = dstB2[i] = 255;
        }
    }

    return 0;
}

}
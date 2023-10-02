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

Matrix blur(Matrix &m, const int radius)
{
    Matrix scratch { PPM::max_dimension };

    double w[Gauss::max_radius] {};
    Gauss::get_weights(radius, w);

    //cache value frequently used, never changed
    const auto dstXsize = m.get_x_size();
    const auto dstYSize = m.get_y_size();

    //pointers for r,g,b in dst matrix
    auto dstR = m.get_R();
    auto dstG = m.get_G();
    auto dstB = m.get_B();

    //non constant pointers so values can be changed
    auto dstR2 = m.get_R_nonconst();
    auto dstG2 = m.get_G_nonconst();
    auto dstB2 = m.get_B_nonconst();

    //pointers for r,g,b scratch matrix
    auto scrR = scratch.get_R();
    auto scrG = scratch.get_G();
    auto scrB = scratch.get_B();

    //non constant pointers so values can be changed
    auto scrR2 = m.get_R_nonconst();
    auto scrG2 = m.get_G_nonconst();
    auto scrB2 = m.get_B_nonconst();

    const auto scrXsize = scratch.get_x_size();

    for (auto x { 0 }; x < dstXsize; x++) {
        for (auto y { 0 }; y < dstYSize; y++) {
            auto r { w[0] * m.r(x, y) }, g { w[0] * m.g(x, y) }, b { w[0] * m.b(x, y) }, n { w[0] };

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
            //scratch.r(x, y) = r / n;
            //scratch.g(x, y) = g / n;
            //scratch.b(x, y) = b / n;
            scrR2[y * scrXsize + x] = r / n;
            scrG2[y * scrXsize + x] = g / n;
            scrB2[y * scrXsize + x] = b / n;
        }
    }

    for (auto x { 0 }; x < dstXsize; x++) {
        for (auto y { 0 }; y < dstYSize; y++) {
            auto r { w[0] * scratch.r(x, y) }, g { w[0] * scratch.g(x, y) }, b { w[0] * scratch.b(x, y) }, n { w[0] };

            for (auto wi { 1 }; wi <= radius; wi++) {
                auto wc { w[wi] };
                auto y2 { y - wi };
                if (y2 >= 0) {
                    //r += wc * scratch.r(x, y2);
                    //g += wc * scratch.g(x, y2);
                    //b += wc * scratch.b(x, y2);
                    r += wc * scrR[y2 * scrXsize + x];
                    g += wc * scrG[y2 * scrXsize + x];
                    b += wc * scrB[y2 * scrXsize + x];
                    n += wc;
                }
                y2 = y + wi;
                if (y2 < dstYSize) {
                    //r += wc * scratch.r(x, y2);
                    //g += wc * scratch.g(x, y2);
                    //b += wc * scratch.b(x, y2);
                    r += wc * scrR[y2 * scrXsize + x];
                    g += wc * scrG[y2 * scrXsize + x];
                    b += wc * scrB[y2 * scrXsize + x];
                    n += wc;
                }
            }
            //m.r(x, y) = r / n;
            //m.g(x, y) = g / n;
            //m.b(x, y) = b / n;

            dstR2[y * dstXsize + x] = r / n;
            dstG2[y * dstXsize + x] = g / n;
            dstB2[y * dstXsize + x] = b / n;
        }
    }

    return 0;
}

Matrix threshold(Matrix &m)
{
    //pointers for r,g,b in dst matrix
    auto dstR = m.get_R();
    auto dstG = m.get_G();
    auto dstB = m.get_B();

    unsigned sum {}, nump { m.get_x_size() * m.get_y_size() };

    for (auto i { 0 }; i < nump; i+=4) {
        //sum += dst.r(i, 0) + dst.g(i, 0) + dst.b(i, 0);
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

}

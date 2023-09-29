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

Matrix blur(Matrix m, const int radius)
{
    Matrix scratch { PPM::max_dimension };
    auto dst { m };

    double w[Gauss::max_radius] {};
    Gauss::get_weights(radius, w);

    //cache value frequently used, never changed
    const auto dstXsize = dst.get_x_size();
    const auto dstYSize = dst.get_y_size();

    std::tuple<unsigned char, unsigned char, unsigned char> rgbTuple;

    for (auto x { 0 }; x < dstXsize; x++) {
        for (auto y { 0 }; y < dstYSize; y++) {
            auto r { w[0] * dst.r(x, y) }, g { w[0] * dst.g(x, y) }, b { w[0] * dst.b(x, y) }, n { w[0] };

            for (auto wi { 1 }; wi <= radius; wi++) {
                auto wc { w[wi] };
                auto x2 { x - wi };
                if (x2 >= 0) {
                    rgbTuple = dst.rgb(x2, y);
                    r += wc * std::get<0>(rgbTuple);
                    g += wc * std::get<1>(rgbTuple);
                    b += wc * std::get<2>(rgbTuple);
                    n += wc;
                }
                x2 = x + wi;
                if (x2 < dstXsize) {
                    rgbTuple = dst.rgb(x2, y);
                    r += wc * std::get<0>(rgbTuple);
                    g += wc * std::get<1>(rgbTuple);
                    b += wc * std::get<2>(rgbTuple);
                    n += wc;
                }
            }
            scratch.r(x, y) = r / n;
            scratch.g(x, y) = g / n;
            scratch.b(x, y) = b / n;
        }
    }

    for (auto x { 0 }; x < dstXsize; x++) {
        for (auto y { 0 }; y < dstYSize; y++) {
            auto r { w[0] * scratch.r(x, y) }, g { w[0] * scratch.g(x, y) }, b { w[0] * scratch.b(x, y) }, n { w[0] };

            for (auto wi { 1 }; wi <= radius; wi++) {
                auto wc { w[wi] };
                auto y2 { y - wi };
                if (y2 >= 0) {
                    rgbTuple = dst.rgb(y2, y);
                    r += wc * std::get<0>(rgbTuple);
                    g += wc * std::get<0>(rgbTuple);
                    b += wc * std::get<0>(rgbTuple);
                    n += wc;
                }
                y2 = y + wi;
                if (y2 < dstYSize) {
                    rgbTuple = dst.rgb(y2, y);
                    r += wc * std::get<0>(rgbTuple);
                    g += wc * std::get<0>(rgbTuple);
                    b += wc * std::get<0>(rgbTuple);
                    n += wc;
                }
            }
            dst.r(x, y) = r / n;
            dst.g(x, y) = g / n;
            dst.b(x, y) = b / n;
        }
    }

    return dst;
}

Matrix threshold(Matrix m)
{
    auto dst { m };
    unsigned sum {}, nump { dst.get_x_size() * dst.get_y_size() };

    for (auto i { 0 }; i < nump; i++) {
        sum += dst.r(i, 0) + dst.g(i, 0) + dst.b(i, 0);
    }

    sum /= nump;

    unsigned psum {};

    for (auto i { 0 }; i < nump; i++) {
        psum = dst.r(i, 0) + dst.g(i, 0) + dst.b(i, 0);
        if (sum > psum) {
            dst.r(i, 0) = dst.g(i, 0) = dst.b(i, 0) = 0;
        } else {
            dst.r(i, 0) = dst.g(i, 0) = dst.b(i, 0) = 255;
        }
    }

    return dst;
}

}

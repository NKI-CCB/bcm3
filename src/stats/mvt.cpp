#include "Utils.h"
#include "checks.h"
#include "mvt.h"

#include <boost/math/special_functions/erf.hpp>
#include <unsupported/Eigen/FFT>

namespace bcm3 {

Real Phi(Real z)
{
    return boost::math::erfc(-z / 1.4142135623731) / 2.0;
}

void CholeskyPermuteReorder(const MatrixReal& sigma, const VectorReal& a, const VectorReal& b, MatrixReal& c, VectorReal& ap, VectorReal& bp)
{
    Real epsilon = -1e-10;

    int n = sigma.cols();
    c = sigma;
    ap = a;
    bp = b;

    VectorReal d = c.diagonal().cwiseMax(VectorReal::Zero(n)).array().sqrt();

    for (int i = 0; i < n; i++) {
        if (d[i] > 0) {
            c.col(i) /= d[i];
            c.row(i) /= d[i];
            ap[i] /= d[i];
            bp[i] /= d[i];
        }
    }

    VectorReal y = VectorReal::Zero(n);
    Real sqtp = sqrt(2.0 * M_PI);
    for (int k = 0; k < n; k++) {
        int im = k;
        Real ckk = 0.0;
        Real dem = 1.0;
        Real s = 0.0;
        Real am = std::numeric_limits<Real>::quiet_NaN();
        Real bm = std::numeric_limits<Real>::quiet_NaN();

        for (int i = k; i < n; i++) {
            if (c(i, i) > std::numeric_limits<Real>::epsilon()) {
                Real cii = sqrt(std::max(c(i, i), 0.0));
                if (i > 0) {
                    // TODO - is this the rigth multiplication?
                    s = c.row(i).segment(0, k - 1) * y.segment(0, k - 1);
                }
                Real ai = (ap[i] - s) / cii;
                Real bi = (bp[i] - s) / cii;
                Real de = Phi(bi) - Phi(ai);
                if (de <= dem) {
                    ckk = cii;
                    dem = de;
                    am = ai;
                    bm = bi;
                    im = i;
                }
            }
        }

        if (im > k) {
            ap[c(im, k)] = ap[c(k, im)];
            bp[c(im, k)] = bp[c(k, im)];

            c(im, im) = c(k, k);

            if (k > 0) {
                VectorReal t = c.row(im).segment(0, k - 1);
                c.row(im).segment(0, k - 1) = c.row(k).segment(0, k - 1);
                c.row(k).segment(0, k - 1) = t;
            }

            if (im < n - 1) {
                VectorReal t = c.col(im).segment(im + 1, n - (im + 1));
                c.col(im).segment(im + 1, n - (im + 1)) = c.col(k).segment(im + 1, n - (im + 1));
                c.col(k).segment(im + 1, n - (im + 1)) = t;
            }

            if (im > k + 1) {
                VectorReal t = c.col(k).segment(k + 1, (im - 1) - (k + 1));
                c.col(k).segment(k + 1, (im - 1) - (k + 1)) = c.row(im).segment(k + 1, (im - 1) - (k + 1)).transpose();
                c.row(im).segment(k + 1, (im - 1) - (k + 1)) = t.transpose();
            }
        }
        if (ckk > epsilon * k) {
            c(k, k) = ckk;
            if (k < n) {
                c.row(k).segment(k + 1, n - (k + 1)).setZero();
                for (int i = k + 1; i < n; i++) {
                    c(i, k) = c(i, k) / ckk;
                    c.row(i).segment(k + 1, i - (k + 1)) -= c(i, k) * c.col(k).segment(k + 1, i - (k + 1)).transpose();
                }
            }
            if (fabs(dem) > epsilon) {
                y[k] = (exp(-(am*am) / 2.0) - exp(-(bm*bm) / 2.0)) / (sqtp * dem);
            } else {
                y[k] = (am + bm) / 2.0;
                if (am < -10.0) {
                    y[k] = bm;
                } else if (bm > 10.0) {
                    y[k] = am;
                }
            }
            c.row(k).segment(0, k) /= ckk;
            ap[k] /= ckk;
            bp[k] /= ckk;
        } else {
            c.col(k).segment(k, n - k).setZero();
            y[k] = (ap[k] + bp[k]) / 2.0;
        }
    }
}

}

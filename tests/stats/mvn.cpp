#define BOOST_TEST_MODULE test_stats_mvn
#include <boost/test/unit_test.hpp>
#include "mvn.h"
#include "ProbabilityDistributions.h"

BOOST_AUTO_TEST_CASE(dmvnorm_1d) {
    VectorReal mu = VectorReal::Zero(1);
    MatrixReal sigma = MatrixReal::Ones(1, 1);
    BOOST_CHECK_CLOSE(bcm3::dmvnormal(VectorReal::Zero(1), mu, sigma), bcm3::PdfNormal(0, mu(0), sigma(0, 0)), 1e-13);
    BOOST_CHECK_CLOSE(bcm3::dmvnormal(VectorReal::Ones(1), mu, sigma), bcm3::PdfNormal(1, mu(0), sigma(0, 0)), 1e-13);
    BOOST_CHECK_CLOSE(bcm3::dmvnormal(-VectorReal::Ones(1), mu, sigma), bcm3::PdfNormal(-1, mu(0), sigma(0, 0)), 1e-13);
}

BOOST_AUTO_TEST_CASE(dmvnorm_2d_identity) {
    VectorReal mu = VectorReal::Zero(2);
    MatrixReal sigma = MatrixReal::Identity(2, 2);
    BOOST_CHECK_CLOSE(bcm3::dmvnormal(VectorReal::Zero(2), mu, sigma), 0.1591549430918953, 1e-12);
    BOOST_CHECK_CLOSE(bcm3::dmvnormal(VectorReal::Ones(2), mu, sigma), 0.05854983152431917, 1e-12);
    BOOST_CHECK_CLOSE(bcm3::dmvnormal(-VectorReal::Ones(2), mu, sigma), 0.05854983152431917, 1e-12);
}

BOOST_AUTO_TEST_CASE(dmvnorm_2d) {
    VectorReal mu = VectorReal::Ones(2) * 1.2e3;
    MatrixReal sigma(2, 2);
    sigma(0, 0) = 234.0;
    sigma(0, 1) = sigma(1, 0) = 42.0;
    sigma(1, 1) = 786;
    BOOST_CHECK_CLOSE(bcm3::dmvnormal(mu, mu, sigma), 0.0003729010642586194, 1e-12);

    VectorReal x(2);
    x(0) = 1000.0;
    x(1) = 1050.0;
    BOOST_CHECK_CLOSE(bcm3::dmvnormal(x, mu, sigma), 6.617956105689106e-45, 1e-12);
}

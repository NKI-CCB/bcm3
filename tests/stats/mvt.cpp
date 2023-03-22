#include <boost/test/unit_test.hpp>
#include "mvt.h"
#include "ProbabilityDistributions.h"

BOOST_AUTO_TEST_CASE(dmvt_1d)
{
	VectorReal mu = VectorReal::Zero(1);
	MatrixReal sigma = MatrixReal::Ones(1, 1);
	BOOST_CHECK_CLOSE(bcm3::dmvt(VectorReal::Zero(1), mu, sigma, 5.0), bcm3::PdfT(0, mu(0), sigma(0, 0), 5.0), 1e-13);
	BOOST_CHECK_CLOSE(bcm3::dmvt(VectorReal::Ones(1), mu, sigma, 5.0), bcm3::PdfT(1, mu(0), sigma(0, 0), 5.0), 1e-13);
	BOOST_CHECK_CLOSE(bcm3::dmvt(-VectorReal::Ones(1), mu, sigma, 5.0), bcm3::PdfT(-1, mu(0), sigma(0, 0), 5.0), 1e-13);
	BOOST_CHECK_CLOSE(bcm3::dmvt(VectorReal::Zero(1), mu, sigma, 5.0, true), bcm3::LogPdfT(0, mu(0), sigma(0, 0), 5.0), 1e-13);
	BOOST_CHECK_CLOSE(bcm3::dmvt(VectorReal::Ones(1), mu, sigma, 5.0, true), bcm3::LogPdfT(1, mu(0), sigma(0, 0), 5.0), 1e-13);
	BOOST_CHECK_CLOSE(bcm3::dmvt(-VectorReal::Ones(1), mu, sigma, 5.0, true), bcm3::LogPdfT(-1, mu(0), sigma(0, 0), 5.0), 1e-13);
}

BOOST_AUTO_TEST_CASE(dmvt_2d_identity)
{
	VectorReal mu = VectorReal::Zero(2);
	MatrixReal sigma = MatrixReal::Identity(2, 2);
	BOOST_CHECK_CLOSE(bcm3::dmvt(VectorReal::Zero(2), mu, sigma, 5.0), 0.1591549430918953, 1e-12);
	BOOST_CHECK_CLOSE(bcm3::dmvt(VectorReal::Ones(2), mu, sigma, 5.0), 0.04901985324897605, 1e-12);
	BOOST_CHECK_CLOSE(bcm3::dmvt(-VectorReal::Ones(2), mu, sigma, 5.0), 0.04901985324897605, 1e-12);
	BOOST_CHECK_CLOSE(bcm3::dmvt(VectorReal::Zero(2), mu, sigma, 5.0, true), -1.837877066409345, 1e-12);
	BOOST_CHECK_CLOSE(bcm3::dmvt(VectorReal::Ones(2), mu, sigma, 5.0, true), -3.015529894583591, 1e-12);
	BOOST_CHECK_CLOSE(bcm3::dmvt(-VectorReal::Ones(2), mu, sigma, 5.0, true), -3.015529894583591, 1e-12);
}

BOOST_AUTO_TEST_CASE(dmvt_2d)
{
	VectorReal mu = VectorReal::Ones(2) * 1.2e3;
	MatrixReal sigma(2, 2);
	sigma(0, 0) = 234.0;
	sigma(0, 1) = sigma(1, 0) = 42.0;
	sigma(1, 1) = 786;
	BOOST_CHECK_CLOSE(bcm3::dmvt(mu, mu, sigma, 5.0), 0.0003729010642586194, 1e-12);
	BOOST_CHECK_CLOSE(bcm3::dmvt(mu, mu, sigma, 5.0, true), -7.894197416778155, 1e-12);

	VectorReal x(2);
	x(0) = 1000.0;
	x(1) = 1050.0;
	BOOST_CHECK_CLOSE(bcm3::dmvt(x, mu, sigma, 5.0), 1.04998038728666e-09, 1e-12);
	BOOST_CHECK_CLOSE(bcm3::dmvt(x, mu, sigma, 5.0, true), -20.67449435172604, 1e-12);
}

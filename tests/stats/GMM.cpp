#include <boost/test/unit_test.hpp>
#include "GMM.h"

BOOST_AUTO_TEST_CASE(GMM_2d)
{
	std::vector<VectorReal> means(2, VectorReal::Zero(2));
	std::vector<MatrixReal> covariances(2, MatrixReal::Identity(2, 2));
	VectorReal weights(2);

	means[0] << -1, -1;
	means[1] << 2, 3;

	covariances[0] << 2, 1, 1, 1;
	covariances[1] << 1, -0.9, -0.9, 1;

	weights << 0.25, 0.75;

	bcm3::GMM gmm;
	gmm.Set(means, covariances, weights);

	VectorReal testpoint(2);
	testpoint << 2, 2;

	BOOST_CHECK_CLOSE(gmm.LogPdf(testpoint), -3.9045912795091535, 1e-12);

	VectorReal responsibilities = gmm.CalculateResponsibilities(testpoint);

	BOOST_CHECK_CLOSE(responsibilities(0), 0.0219370092578219, 1e-12);
	BOOST_CHECK_CLOSE(responsibilities(1), 0.9780629907421782, 1e-12);
	
}

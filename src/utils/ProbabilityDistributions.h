#pragma once

namespace bcm3 {
	
Real PdfBeta(Real x, Real a, Real b);
Real PdfBetaPrime(Real x, Real a, Real b, Real scale);
Real PdfExponential(Real x, Real lambda);
Real PdfGamma(Real x, Real shape, Real scale);
Real PdfNormal(Real x, Real mu, Real sigma);
Real PdfT(Real x, Real mu, Real sigma, Real nu);
Real PdfT(Real x, Real mu, Real sigma, Real nu, Real C);
Real PdfGPD(Real x, Real u, Real xi, Real beta);

Real LogPdfBetaPrecision(Real x, Real mu, Real phi);
Real LogPdfBetaPrime(Real x, Real a, Real b, Real scale);
Real LogPdfExponential(Real x, Real lambda);
Real LogPdfNormal(Real x, Real mu, Real sigma, bool skip_na = false);
Real LogPdfTruncatedNormal(Real x, Real mu, Real sigma, Real a, Real b, bool skip_na = false);
Real LogPdfT(Real x, Real mu, Real sigma, Real nu, bool skip_na = false);
Real LogPdfT(Real x, Real mu, Real sigma, Real nu, Real logC, bool skip_na = false);
Real LogPdfT_CalcC(Real nu);
Real LogPdfTnu3(Real x, Real mu, Real sigma, bool skip_na = false);
Real LogTruncatedPdfTnu3(Real x, Real mu, Real sigma, Real a, Real b, bool skip_na = false);
Real LogPdfGPD(Real x, Real u, Real xi, Real beta);

Real CdfBeta(Real x, Real a, Real b);
Real CdfCauchy(Real x, Real scale);
Real CdfExponential(Real x, Real lambda);
Real CdfGamma(Real x, Real shape, Real scale);
Real CdfHalfCauchy(Real x, Real scale);
Real CdfNormal(Real x, Real mu, Real sigma);
Real CdfT(Real x, Real mu, Real sigma, Real nu);
Real CdfGPD(Real x, Real u, Real xi, Real beta);
Real CdfTruncatedNormal(Real x, Real mu, Real sigma, Real a, Real b);

Real QuantileBeta(Real p, Real a, Real b);
Real QuantileExponential(Real p, Real lambda);
Real QuantileGamma(Real p, Real shape, Real scale);
Real QuantileHalfCauchy(Real p, Real scale);
Real QuantileNormal(Real p, Real mu, Real sigma);
Real QuantileT(Real p, Real mu, Real sigma, Real nu);
Real QuantileUniform(Real p, Real a, Real b);
Real QuantileGPD(Real x, Real u, Real xi, Real beta);

}

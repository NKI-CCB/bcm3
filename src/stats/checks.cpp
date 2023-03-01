#include "Utils.h"
#include "checks.h"

namespace bcm3
{

bool is_positive_semi_definite(const MatrixReal& A)
{
    if (!A.isApprox(A.transpose())) {
        return false;
    }
    const auto ldlt = A.template selfadjointView<Eigen::Upper>().ldlt();
    if (ldlt.info() == Eigen::NumericalIssue || !ldlt.isPositive()) {
        return false;
    }
    return true;
}

}

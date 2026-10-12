#pragma once

namespace bcm3 {

// Returns match(i) = the column assigned to row i.
// Returns nullopt if not every row can be assigned.
std::optional<Eigen::VectorXi> RectangularHungarian(const Eigen::MatrixXd& cost);

}

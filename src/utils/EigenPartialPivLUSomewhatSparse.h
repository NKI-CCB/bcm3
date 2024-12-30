#pragma once

template<typename MatrixType>
class PartialPivLUExtended : public Eigen::PartialPivLU< MatrixType >
{
	typedef typename MatrixType::Index Index;
	typedef typename MatrixType::RealScalar RealScalar;
	typedef typename MatrixType::Scalar Scalar;

public:
	template<typename InputType>
	PartialPivLUExtended& compute_select(const Eigen::EigenBase<InputType>& matrix)
	{
		if (matrix.rows() == 3) {
			ASSERT(matrix.cols() == 3);
			inverse.resize(3, 3);
			Eigen::internal::compute_inverse<MatrixType, MatrixType, 3>::run(matrix, inverse);
			return *this;
		} else {
			return compute_optimized(matrix);
		}
	}

	template<typename RhsType, typename DstType>
	void apply_select(const RhsType& b, DstType& x)
	{
		if (b.size() == 3) {
			x.noalias() = inverse * b;
		} else {
			x.noalias() = Eigen::PartialPivLU< MatrixType >::solve(b);
		}
	}

	template<typename InputType>
	PartialPivLUExtended& compute_optimized(const Eigen::EigenBase<InputType>& matrix)
	{
		this->m_lu = matrix.derived();

		this->check_template_parameters();

		if (this->m_lu.cols() > 0)
			this->m_l1_norm = this->m_lu.cwiseAbs().colwise().sum().maxCoeff();
		else
			this->m_l1_norm = RealScalar(0);

		eigen_assert(this->m_lu.rows() == this->m_lu.cols() && "PartialPivLU is only for square (and moreover invertible) matrices");
		const Index size = this->m_lu.rows();

		this->m_rowsTranspositions.resize(size);
		typename Eigen::PartialPivLU<MatrixType>::TranspositionType::StorageIndex nb_transpositions;

		// Following code block is copied and modified from unblocked_lu(MatrixTypeRef & lu, PivIndex * row_transpositions, PivIndex & nb_transpositions)
		{
			typedef Eigen::internal::scalar_score_coeff_op<Scalar> Scoring;
			typedef typename Scoring::result_type Score;
			const Index rows = this->m_lu.rows();
			const Index cols = this->m_lu.cols();
			const Index size = (std::min)(rows, cols);
			static const int RRows = MatrixType::SizeAtCompileTime == 2 ? 1 : Eigen::Dynamic;
			static const int RCols = MatrixType::SizeAtCompileTime == 2 ? 1 : Eigen::Dynamic;

			nb_transpositions = 0;
			for (Index k = 0; k < size; ++k) {
				int rrows = Eigen::internal::convert_index<int>(rows - k - 1);
				int rcols = Eigen::internal::convert_index<int>(cols - k - 1);

				Index row_of_biggest_in_col;
				Score biggest_in_corner = this->m_lu.col(k).tail(rows - k).unaryExpr(Scoring()).maxCoeff(&row_of_biggest_in_col);
				row_of_biggest_in_col += k;
				this->m_rowsTranspositions[k] = row_of_biggest_in_col;

				if (biggest_in_corner != Score(0)) {
					if (k != row_of_biggest_in_col) {
						this->m_lu.row(k).swap(this->m_lu.row(row_of_biggest_in_col));
						++nb_transpositions;
					}

					Scalar inv_coeff = 1.0 / this->m_lu.coeff(k, k);
					this->m_lu.col(k).tail(Eigen::fix<RRows>(rrows)) *= inv_coeff;
				}

#if 0
				if (k < rows - 1) {
					m_lu.bottomRightCorner(Eigen::fix<RRows>(rrows), Eigen::fix<RCols>(rcols)).noalias() -= m_lu.col(k).tail(Eigen::fix<RRows>(rrows)) * m_lu.row(k).tail(Eigen::fix<RCols>(rcols));
				}
#else
				// This is the part that's modified - when matrices have quite a few 0's, we can skip many of the entries
				for (Index j = k + 1; j < cols; j++) {
					Scalar a_kj = this->m_lu.coeff(k, j);
					if (a_kj != 0.0) {
						this->m_lu.col(j).tail(Eigen::fix<RRows>(rrows)).noalias() -= a_kj * this->m_lu.col(k).tail(Eigen::fix<RRows>(rrows));
					}
				}
#endif
			}
		}

		this->m_det_p = (nb_transpositions % 2) ? -1 : 1;
		this->m_p = this->m_rowsTranspositions;
		this->m_isInitialized = true;
		return *this;
	}

	MatrixType inverse;
};

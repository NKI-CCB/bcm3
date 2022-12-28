#pragma once

template<typename MatrixType>
class PartialPivLUExtended : public Eigen::PartialPivLU< MatrixType >
{
	typedef Eigen::Ref<MatrixType> MatrixTypeRef;

public:
#if 0
	static Index unblocked_lu(MatrixTypeRef& lu, PivIndex* row_transpositions, PivIndex& nb_transpositions)
	{
		typedef scalar_score_coeff_op<Scalar> Scoring;
		typedef typename Scoring::result_type Score;
		const Index rows = lu.rows();
		const Index cols = lu.cols();
		const Index size = (std::min)(rows, cols);
		// For small compile-time matrices it is worth processing the last row separately:
		//  speedup: +100% for 2x2, +10% for others.
		const Index endk = UnBlockedAtCompileTime ? size - 1 : size;
		nb_transpositions = 0;
		Index first_zero_pivot = -1;
		for (Index k = 0; k < endk; ++k)
		{
			int rrows = internal::convert_index<int>(rows - k - 1);
			int rcols = internal::convert_index<int>(cols - k - 1);

			Index row_of_biggest_in_col;
			Score biggest_in_corner
				= lu.col(k).tail(rows - k).unaryExpr(Scoring()).maxCoeff(&row_of_biggest_in_col);
			row_of_biggest_in_col += k;

			row_transpositions[k] = PivIndex(row_of_biggest_in_col);

			if (biggest_in_corner != Score(0))
			{
				if (k != row_of_biggest_in_col)
				{
					lu.row(k).swap(lu.row(row_of_biggest_in_col));
					++nb_transpositions;
				}

				lu.col(k).tail(fix<RRows>(rrows)) /= lu.coeff(k, k);
			} else if (first_zero_pivot == -1)
			{
				// the pivot is exactly zero, we record the index of the first pivot which is exactly 0,
				// and continue the factorization such we still have A = PLU
				first_zero_pivot = k;
			}

			if (k < rows - 1)
				lu.bottomRightCorner(fix<RRows>(rrows), fix<RCols>(rcols)).noalias() -= lu.col(k).tail(fix<RRows>(rrows)) * lu.row(k).tail(fix<RCols>(rcols));
		}

		// special handling of the last entry
		if (UnBlockedAtCompileTime)
		{
			Index k = endk;
			row_transpositions[k] = PivIndex(k);
			if (Scoring()(lu(k, k)) == Score(0) && first_zero_pivot == -1)
				first_zero_pivot = k;
		}

		return first_zero_pivot;
	}

	template<typename InputType>
	PartialPivLUExtended& compute_optimized_blocked(const EigenBase<InputType>& matrix)
	{
		m_lu = matrix.derived();

		check_template_parameters();

		if (m_lu.cols() > 0)
			m_l1_norm = m_lu.cwiseAbs().colwise().sum().maxCoeff();
		else
			m_l1_norm = RealScalar(0);

		eigen_assert(m_lu.rows() == m_lu.cols() && "PartialPivLU is only for square (and moreover invertible) matrices");
		const Index size = m_lu.rows();

		m_rowsTranspositions.resize(size);
		typename TranspositionType::StorageIndex nb_transpositions;

		// Following code block is copied and modified from blocked_lu(Index rows, Index cols, Scalar* lu_data, Index luStride, PivIndex* row_transpositions, PivIndex& nb_transpositions, Index maxBlockSize=256)
		{
			MatrixTypeRef lu = MatrixType::Map(lu_data, rows, cols, OuterStride<>(luStride));

			const Index size = (std::min)(rows, cols);

			// if the matrix is too small, no blocking:
			if (UnBlockedAtCompileTime || size <= UnBlockedBound)
			{
				return unblocked_lu(lu, row_transpositions, nb_transpositions);
			}

			// automatically adjust the number of subdivisions to the size
			// of the matrix so that there is enough sub blocks:
			Index blockSize;
			{
				blockSize = size / 8;
				blockSize = (blockSize / 16) * 16;
				blockSize = (std::min)((std::max)(blockSize, Index(8)), maxBlockSize);
			}

			nb_transpositions = 0;
			Index first_zero_pivot = -1;
			for (Index k = 0; k < size; k += blockSize)
			{
				Index bs = (std::min)(size - k, blockSize); // actual size of the block
				Index trows = rows - k - bs; // trailing rows
				Index tsize = size - k - bs; // trailing size

				// partition the matrix:
				//                          A00 | A01 | A02
				// lu  = A_0 | A_1 | A_2 =  A10 | A11 | A12
				//                          A20 | A21 | A22
				BlockType A_0 = lu.block(0, 0, rows, k);
				BlockType A_2 = lu.block(0, k + bs, rows, tsize);
				BlockType A11 = lu.block(k, k, bs, bs);
				BlockType A12 = lu.block(k, k + bs, bs, tsize);
				BlockType A21 = lu.block(k + bs, k, trows, bs);
				BlockType A22 = lu.block(k + bs, k + bs, trows, tsize);

				PivIndex nb_transpositions_in_panel;
				// recursively call the blocked LU algorithm on [A11^T A21^T]^T
				// with a very small blocking size:
				Index ret = blocked_lu(trows + bs, bs, &lu.coeffRef(k, k), luStride,
					row_transpositions + k, nb_transpositions_in_panel, 16);
				if (ret >= 0 && first_zero_pivot == -1)
					first_zero_pivot = k + ret;

				nb_transpositions += nb_transpositions_in_panel;
				// update permutations and apply them to A_0
				for (Index i = k; i < k + bs; ++i)
				{
					Index piv = (row_transpositions[i] += internal::convert_index<PivIndex>(k));
					A_0.row(i).swap(A_0.row(piv));
				}

				if (trows)
				{
					// apply permutations to A_2
					for (Index i = k; i < k + bs; ++i)
						A_2.row(i).swap(A_2.row(row_transpositions[i]));

					// A12 = A11^-1 A12
					A11.template triangularView<UnitLower>().solveInPlace(A12);

					A22.noalias() -= A21 * A12;
				}
			}

			typedef Eigen::internal::scalar_score_coeff_op<Scalar> Scoring;
			typedef typename Scoring::result_type Score;
			const Index rows = m_lu.rows();
			const Index cols = m_lu.cols();
			const Index size = (std::min)(rows, cols);
			static const int RRows = SizeAtCompileTime == 2 ? 1 : Eigen::Dynamic;
			static const int RCols = SizeAtCompileTime == 2 ? 1 : Eigen::Dynamic;

			nb_transpositions = 0;
			for (Index k = 0; k < size; ++k) {
				int rrows = Eigen::internal::convert_index<int>(rows - k - 1);
				int rcols = Eigen::internal::convert_index<int>(cols - k - 1);

				Index row_of_biggest_in_col;
				Score biggest_in_corner = m_lu.col(k).tail(rows - k).unaryExpr(Scoring()).maxCoeff(&row_of_biggest_in_col);
				row_of_biggest_in_col += k;
				m_rowsTranspositions[k] = row_of_biggest_in_col;

				if (biggest_in_corner != Score(0)) {
					if (k != row_of_biggest_in_col) {
						m_lu.row(k).swap(m_lu.row(row_of_biggest_in_col));
						++nb_transpositions;
					}

					Scalar inv_coeff = 1.0 / m_lu.coeff(k, k);
					m_lu.col(k).tail(Eigen::fix<RRows>(rrows)) *= inv_coeff;
				}

#if 0
				if (k < rows - 1) {
					m_lu.bottomRightCorner(Eigen::fix<RRows>(rrows), Eigen::fix<RCols>(rcols)).noalias() -= m_lu.col(k).tail(Eigen::fix<RRows>(rrows)) * m_lu.row(k).tail(Eigen::fix<RCols>(rcols));
				}
#else
				// This is the part that's modified - when matrices have quite a few 0's, we can skip many of the entries
				for (Index j = k + 1; j < cols; j++) {
					Scalar a_kj = m_lu.coeff(k, j);
					if (a_kj != 0.0) {
						m_lu.col(j).tail(Eigen::fix<RRows>(rrows)).noalias() -= a_kj * m_lu.col(k).tail(Eigen::fix<RRows>(rrows));
					}
				}
#endif
			}
		}

		m_det_p = (nb_transpositions % 2) ? -1 : 1;
		m_p = m_rowsTranspositions;
		m_isInitialized = true;
		return *this;
	}
#endif
	template<typename InputType>
	PartialPivLUExtended& compute_optimized(const EigenBase<InputType>& matrix)
	{
		m_lu = matrix.derived();

		check_template_parameters();

		if (m_lu.cols() > 0)
			m_l1_norm = m_lu.cwiseAbs().colwise().sum().maxCoeff();
		else
			m_l1_norm = RealScalar(0);

		eigen_assert(m_lu.rows() == m_lu.cols() && "PartialPivLU is only for square (and moreover invertible) matrices");
		const Index size = m_lu.rows();

		m_rowsTranspositions.resize(size);
		typename TranspositionType::StorageIndex nb_transpositions;

		// Following code block is copied and modified from unblocked_lu(MatrixTypeRef & lu, PivIndex * row_transpositions, PivIndex & nb_transpositions)
		{
			typedef Eigen::internal::scalar_score_coeff_op<Scalar> Scoring;
			typedef typename Scoring::result_type Score;
			const Index rows = m_lu.rows();
			const Index cols = m_lu.cols();
			const Index size = (std::min)(rows, cols);
			static const int RRows = SizeAtCompileTime == 2 ? 1 : Eigen::Dynamic;
			static const int RCols = SizeAtCompileTime == 2 ? 1 : Eigen::Dynamic;

			nb_transpositions = 0;
			for (Index k = 0; k < size; ++k) {
				int rrows = Eigen::internal::convert_index<int>(rows - k - 1);
				int rcols = Eigen::internal::convert_index<int>(cols - k - 1);

				Index row_of_biggest_in_col;
				Score biggest_in_corner = m_lu.col(k).tail(rows - k).unaryExpr(Scoring()).maxCoeff(&row_of_biggest_in_col);
				row_of_biggest_in_col += k;
				m_rowsTranspositions[k] = row_of_biggest_in_col;

				if (biggest_in_corner != Score(0)) {
					if (k != row_of_biggest_in_col) {
						m_lu.row(k).swap(m_lu.row(row_of_biggest_in_col));
						++nb_transpositions;
					}

					Scalar inv_coeff = 1.0 / m_lu.coeff(k, k);
					m_lu.col(k).tail(Eigen::fix<RRows>(rrows)) *= inv_coeff;
				}

#if 0
				if (k < rows - 1) {
					m_lu.bottomRightCorner(Eigen::fix<RRows>(rrows), Eigen::fix<RCols>(rcols)).noalias() -= m_lu.col(k).tail(Eigen::fix<RRows>(rrows)) * m_lu.row(k).tail(Eigen::fix<RCols>(rcols));
				}
#else
				// This is the part that's modified - when matrices have quite a few 0's, we can skip many of the entries
				for (Index j = k + 1; j < cols; j++) {
					Scalar a_kj = m_lu.coeff(k, j);
					if (a_kj != 0.0) {
						m_lu.col(j).tail(Eigen::fix<RRows>(rrows)).noalias() -= a_kj * m_lu.col(k).tail(Eigen::fix<RRows>(rrows));
					}
				}
#endif
			}
		}

		m_det_p = (nb_transpositions % 2) ? -1 : 1;
		m_p = m_rowsTranspositions;
		m_isInitialized = true;
		return *this;
	}
};

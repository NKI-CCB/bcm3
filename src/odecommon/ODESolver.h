#pragma once

#if ODE_SINGLE_PRECISION
	typedef float OdeReal;
	typedef Eigen::VectorXf OdeVectorReal;
	typedef Eigen::MatrixXf OdeMatrixReal;
#else
	typedef double OdeReal;
	typedef Eigen::VectorXd OdeVectorReal;
	typedef Eigen::MatrixXd OdeMatrixReal;
#endif

class ODESolver
{
public:
	typedef std::function<bool(OdeReal, const OdeReal*, OdeReal*, void*)> TDeriviativeFunction;
	typedef std::function<bool(OdeReal, const OdeReal*, const OdeReal*, OdeReal**, void*)> TJacobianFunction;
	typedef std::function<Real(OdeReal, void*)> TDiscontinuityCallback;

	ODESolver();
	virtual ~ODESolver();

	virtual bool Initialize(size_t N, void* user);
	bool SetTolerance(OdeReal relative, OdeReal absolute);
	void SetUserData(void* user);
	void SetDiscontinuity(OdeReal time, TDiscontinuityCallback cb, void* user);
	void ResetDiscontinuity();

	void SetDerivativeFunction(TDeriviativeFunction f);
	void SetJacobianFunction(TJacobianFunction f);

	virtual bool Simulate(const Real* initial_conditions, const VectorReal& timepoints, MatrixReal& output, bool verbose = false) = 0;
	virtual OdeReal get_y(size_t i) { return std::numeric_limits<OdeReal>::quiet_NaN(); }
	virtual void set_y(size_t i, OdeReal y) {}

protected:
	size_t N;
	void* user_data;
	Real discontinuity_time;
	TDiscontinuityCallback discontinuity_cb;
	void* discontinuity_user;
	TDeriviativeFunction derivative;
	TJacobianFunction jacobian;

	Real abstol;
	Real reltol;
};

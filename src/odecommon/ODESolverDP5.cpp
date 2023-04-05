#include "Utils.h"
#include "ODESolverDP5.h"
#include <new>
#include <xmmintrin.h>

#define USE_SSE2 1

static const OdeReal MIN_DT = 1e-3;
static const int ATTEMPTS = 5;
static const OdeReal MIN_SCALE_FACTOR = 0.2;
static const OdeReal MAX_SCALE_FACTOR = 5.0;
static const unsigned int MAX_STEPS = 20000;

ODESolverDP5::ODESolverDP5()
	: ytmp(nullptr)
	, yn(nullptr)
{
	for (int i = 0; i < 7; i++) {
		k[i] = nullptr;
	}
}

ODESolverDP5::~ODESolverDP5()
{
	if (ytmp != nullptr) {
		for (int i = 0; i < 7; i++) {
			delete[] k[i];
		}
		delete[] ytmp;
		delete[] yn;
	}
}

bool ODESolverDP5::Initialize(size_t N, void* user)
{
	if (!ODESolver::Initialize(N, user)) {
		return false;
	}

	if (ytmp != nullptr) {
		for (int i = 0; i < 7; i++) {
			delete[] k[i];
		}
		delete[] ytmp;
		delete[] yn;
	}

	try {
		size_t allocsize = N;
		if (N % 2 == 1) {
			allocsize++;
		}
		for (int i = 0; i < 7; i++) {
			k[i] = new OdeReal[allocsize];
			memset(k[i], 0, allocsize * sizeof(OdeReal));
		}
		ytmp = new OdeReal[allocsize];
		yn = new OdeReal[allocsize];
		memset(ytmp, 0, allocsize * sizeof(OdeReal));
		memset(yn, 0, allocsize * sizeof(OdeReal));
	} catch (std::bad_alloc& ba) {
		LOGERROR("Memory allocation failure: %s", ba.what());
		return false;
	}

	return true;
}

bool ODESolverDP5::Simulate(const OdeReal* initial_conditions, const OdeVectorReal& timepoints, OdeMatrixReal& output, bool verbose)
{
	ASSERT(ytmp != nullptr);
	
	if (timepoints.size() == 0) {
		LOGERROR("No timepoints requested");
		return false;
	}

	output.setConstant(N, timepoints.size(), std::numeric_limits<OdeReal>::quiet_NaN());

	size_t ti = 0;
	if (timepoints(0) < std::numeric_limits<OdeReal>::epsilon()) {
		for (size_t i = 0; i < N; i++) {
			output(i, 0) = initial_conditions[i];
		}
		ti++;
	}

	OdeReal t = 0.0;
	OdeReal dt = 0.1;

	for (size_t i = 0; i < N; i++) {
		yn[i] = initial_conditions[i];
	}

	derivative(t, yn, k[0], user_data);

	unsigned int steps = 0;
	
	while (1) {
		// Temporarily decrease timestep if we're approaching a discontinuity.
		OdeReal cur_dt;
		OdeReal next_dt = dt;
		bool approaching_discontinuity = false;
		if (discontinuity_time == discontinuity_time) {
			if (discontinuity_time < t + dt) {
				cur_dt = discontinuity_time - t;
				approaching_discontinuity = true;
			} else {
				cur_dt = dt;
			}
		}
		
		bool succeeded = false;
		for (int i = 0; i < ATTEMPTS; i++) {
			OdeReal maxdiff = ApplyRK(t, cur_dt);
			if (std::isnan(maxdiff) || maxdiff == -std::numeric_limits<OdeReal>::infinity()) {
				// LOGERROR(ODESolverDopri, IntegrateImpl, "NaN in error calculation");
				return false;
			}

			if (verbose) {
				LOG("t=%g - dt=%g maxdiff=%g dctime=%g%s", t, cur_dt, maxdiff, discontinuity_time, approaching_discontinuity ? " approaching dc" : "");
			}

			// From Hairer I; section II.4, p167
			if (maxdiff > 1.1) {
				if (cur_dt == MIN_DT) {
					break;
				} else {
					OdeReal scale = (OdeReal)0.9 * pow(maxdiff, (OdeReal)-0.2);
					scale = std::max(MIN_SCALE_FACTOR, scale);
					cur_dt *= scale;
					if (cur_dt < MIN_DT) {
						cur_dt = MIN_DT;
					}
					if (discontinuity_time == discontinuity_time) {
						if (approaching_discontinuity && t + cur_dt < discontinuity_time) {
							// Time step decreased enough so that we're no longer approaching the discontinuity in this time step
							approaching_discontinuity = false;
						}
					}
				}
			} else if (maxdiff < (OdeReal)0.5) {
				if (!approaching_discontinuity) {
					maxdiff = std::max(maxdiff, (OdeReal)1e-5);
					OdeReal scale = (OdeReal)0.9 * pow(maxdiff, (OdeReal)-0.2);
					scale = std::min(MAX_SCALE_FACTOR, scale);
					next_dt = cur_dt * scale;
				}
				succeeded = true;
				break;
			} else {
				next_dt = cur_dt;
				succeeded = true;
				break;
			}
		}

		if (!succeeded) {
			//LOGERROR(ODESolverDopri, IntegrateImpl, "Time step adaptation did not converge in %d steps", ATTEMPTS);
			return false;
		}

		// Interpolate any timepoints we may have passed
		OdeReal target_t = t + cur_dt;
		while (target_t >= timepoints(ti)) {
			OdeReal theta = timepoints(ti) - t;
			if (theta >= 1.0) { // theta should not be bigger than 1.0; if it is than it should be only a rounding error
				for (size_t i = 0; i < N; i++) {
					output(i, ti) = ytmp[i];
				}
			} else {
				// From Hairer I; section II.5, p179
				OdeReal thetaSq = theta * theta;
				OdeReal b1 = theta * ((OdeReal)1.0 + theta * ((OdeReal)-2.7854166666666669 + theta * ((OdeReal)2.8861111111111111 + theta * ((OdeReal)-1.0095486111111112))));
				OdeReal b3 = (OdeReal)33.33333333333333 * thetaSq * (0.11363881401617251 + theta * ((OdeReal)-0.1682659478885894 + theta * (OdeReal)0.068104222821203958));
				OdeReal b4 = (OdeReal)-2.5 * thetaSq * ((OdeReal)0.675 + theta * ((OdeReal)-1.8 + theta * ((OdeReal)0.8645833333333333)));
				OdeReal b5 = (OdeReal)21.491745283018869 * thetaSq * (-0.012 + theta * ((OdeReal)0.058666666666666666 + theta * ((OdeReal)-0.06166666666666666)));
				OdeReal b6 = (OdeReal)-3.1428571428571428 * thetaSq * (-0.3 + theta * ((OdeReal)0.9666666666666666 + theta * ((OdeReal)-0.7083333333333333)));

				for (size_t i = 0; i < N; i++) {
					output(i, ti) = yn[i] + cur_dt * (b1 * k[0][i] +
													  b3 * k[2][i] +
													  b4 * k[3][i] +
													  b5 * k[4][i] +
													  b6 * k[5][i]);
				}
			}
			ti++;
			if (ti == timepoints.size()) {
				break;
			}
		}
		if (ti == timepoints.size()) {
			// All done; no need to finish the step
			break;
		}

		// Advance solution
		memcpy(k[0], k[6], N * sizeof(OdeReal));
		memcpy(yn, ytmp, N * sizeof(OdeReal));

		if (approaching_discontinuity) {
			// Might as well set the time exactly, although t + cur_dt should equal discontinuity_time
			t = discontinuity_time;

			// Call the callback and reset the system if necessary
			discontinuity_time = discontinuity_cb(t, discontinuity_user);
			derivative(t, yn, k[0], user_data);
			next_dt = 0.1;

			steps = 0;
		} else {
			t += cur_dt;
		}

		dt = next_dt;
		steps++;
		if (steps == MAX_STEPS) {
			//LOGERROR(ODESolverDopri, IntegrateImpl, "MAXSTEPS reached, bailing out");
			return false;
		}
	}

	return true;
}

OdeReal ODESolverDP5::get_y(size_t i)
{
	return yn[i];
}

void ODESolverDP5::set_y(size_t i, OdeReal y)
{
	yn[i] = y;
}

OdeReal ODESolverDP5::ApplyRK(OdeReal t, OdeReal cur_dt)
{
#if USE_SSE2
	__m128d* k1p = (__m128d*)k[0];
	__m128d* k2p = (__m128d*)k[1];
	__m128d* k3p = (__m128d*)k[2];
	__m128d* k4p = (__m128d*)k[3];
	__m128d* k5p = (__m128d*)k[4];
	__m128d* k6p = (__m128d*)k[5];
	__m128d* k7p = (__m128d*)k[6];
	__m128d* ytmpp = (__m128d*)ytmp;
	__m128d* ynp = (__m128d*)yn;
	__m128d cur_dt_sse = _mm_set1_pd(cur_dt);
	size_t iters = (N+1)>>1;
#endif

	// K2
#if USE_SSE2
	__m128d f = _mm_mul_pd(cur_dt_sse, _mm_set1_pd(0.2));
	for (size_t i = 0; i < iters; i++) {
		ytmpp[i] = _mm_add_pd(ynp[i], _mm_mul_pd(f, k1p[i]));
	}
#else
	for (size_t i = 0; i < N; i++) {
		ytmp[i] = yn[i] + cur_dt * (OdeReal)0.2 * k1[i];
	}
#endif
	derivative(t + (OdeReal)0.2 * cur_dt, ytmp, k[1], user_data);

	// K3
#if USE_SSE2
	for (size_t i = 0; i < iters; i++) {
		__m128d t;
		t =					_mm_mul_pd(_mm_set1_pd(0.075), k1p[i]);
		t = _mm_add_pd(t,	_mm_mul_pd(_mm_set1_pd(0.225), k2p[i]));
		ytmpp[i] = _mm_add_pd(ynp[i], _mm_mul_pd(cur_dt_sse, t));
	}
#else
	for (size_t i = 0; i < N; i++) {
		ytmp[i] = yn[i] + cur_dt * (
			+ (OdeReal)0.075 * k1[i]
			+ (OdeReal)0.225 * k2[i]);
	}
#endif
	derivative(t + (OdeReal)0.3 * cur_dt, ytmp, k[2], user_data);

	// K4
#if USE_SSE2
	for (size_t i = 0; i < iters; i++) {
		__m128d t;
		t =					_mm_mul_pd(_mm_set1_pd( 0.97777777777777777777777777777778), k1p[i]);
		t = _mm_add_pd(t,	_mm_mul_pd(_mm_set1_pd(-3.7333333333333333333333333333333 ), k2p[i]));
		t = _mm_add_pd(t,	_mm_mul_pd(_mm_set1_pd( 3.5555555555555555555555555555556 ), k3p[i]));
		ytmpp[i] = _mm_add_pd(ynp[i], _mm_mul_pd(cur_dt_sse, t));
	}
#else
	for (size_t i = 0; i < N; i++) {
		ytmp[i] = yn[i] + cur_dt * (
			+ (OdeReal)0.97777777777777777777777777777778 * k1[i]
			- (OdeReal)3.7333333333333333333333333333333  * k2[i]
			+ (OdeReal)3.5555555555555555555555555555556  * k3[i]);
	}
#endif
	derivative(t + (OdeReal)0.8 * cur_dt, ytmp, k[3], user_data);

	// K5
#if USE_SSE2
	for (size_t i = 0; i < iters; i++) {
		__m128d t;
		t =					_mm_mul_pd(_mm_set1_pd(  2.9525986892242036274958085657674 ), k1p[i]);
		t = _mm_add_pd(t,	_mm_mul_pd(_mm_set1_pd(-11.595793324188385916780978509374  ), k2p[i]));
		t = _mm_add_pd(t,	_mm_mul_pd(_mm_set1_pd(  9.8228928516994360615759792714525 ), k3p[i]));
		t = _mm_add_pd(t,	_mm_mul_pd(_mm_set1_pd( -0.29080932784636488340192043895748), k4p[i]));
		ytmpp[i] = _mm_add_pd(ynp[i], _mm_mul_pd(cur_dt_sse, t));
	}
#else
	for (size_t i = 0; i < N; i++) {
		ytmp[i] = yn[i] + cur_dt * (
			+ (OdeReal)2.9525986892242036274958085657674  * k1[i]
			- (OdeReal)11.595793324188385916780978509374  * k2[i]
			+ (OdeReal)9.8228928516994360615759792714525  * k3[i]
			- (OdeReal)0.29080932784636488340192043895748 * k4[i]);
	}
#endif
	derivative(t + (OdeReal)0.88888888888888888888888888888889 * cur_dt, ytmp, k[4], user_data);

	// K6
#if USE_SSE2
	for (size_t i = 0; i < iters; i++) {
		__m128d t;
		t =					_mm_mul_pd(_mm_set1_pd(  2.8462752525252525252525252525253 ), k1p[i]);
		t = _mm_add_pd(t,	_mm_mul_pd(_mm_set1_pd(-10.757575757575757575757575757576  ), k2p[i]));
		t = _mm_add_pd(t,	_mm_mul_pd(_mm_set1_pd(  8.9064227177434724604535925290642 ), k3p[i]));
		t = _mm_add_pd(t,	_mm_mul_pd(_mm_set1_pd(  0.27840909090909090909090909090909), k4p[i]));
		t = _mm_add_pd(t,	_mm_mul_pd(_mm_set1_pd( -0.27353130360205831903945111492281), k5p[i]));
		ytmpp[i] = _mm_add_pd(ynp[i], _mm_mul_pd(cur_dt_sse, t));
	}
#else
	for (size_t i = 0; i < N; i++) {
		ytmp[i] = yn[i] + cur_dt * (
			+ (OdeReal)2.8462752525252525252525252525253  * k1[i]
			- (OdeReal)10.757575757575757575757575757576  * k2[i]
			+ (OdeReal)8.9064227177434724604535925290642  * k3[i]
			+ (OdeReal)0.27840909090909090909090909090909 * k4[i]
			- (OdeReal)0.27353130360205831903945111492281 * k5[i]);
	}
#endif
	derivative(t + cur_dt, ytmp, k[5], user_data);

	// 5th order accurate
#if USE_SSE2
	for (size_t i = 0; i < iters; i++) {
		__m128d t;
		t =					_mm_mul_pd(_mm_set1_pd( 0.09114583333333333333333333333333), k1p[i]);
		t = _mm_add_pd(t,	_mm_mul_pd(_mm_set1_pd( 0.44923629829290206648697214734951), k3p[i]));
		t = _mm_add_pd(t,	_mm_mul_pd(_mm_set1_pd( 0.65104166666666666666666666666667), k4p[i]));
		t = _mm_add_pd(t,	_mm_mul_pd(_mm_set1_pd(-0.32237617924528301886792452830189), k5p[i]));
		t = _mm_add_pd(t,	_mm_mul_pd(_mm_set1_pd( 0.13095238095238095238095238095238), k6p[i]));
		ytmpp[i] = _mm_add_pd(ynp[i], _mm_mul_pd(cur_dt_sse, t));
	}
#else
	for (size_t i = 0; i < N; i++) {
		ytmp[i] = yn[i] + cur_dt * (
			+ (OdeReal)0.09114583333333333333333333333333 * k1[i]
			+ (OdeReal)0.44923629829290206648697214734951 * k3[i]
			+ (OdeReal)0.65104166666666666666666666666667 * k4[i]
			- (OdeReal)0.32237617924528301886792452830189 * k5[i]
			+ (OdeReal)0.13095238095238095238095238095238 * k6[i]);
	}
#endif

	// K7
	derivative(t + cur_dt, ytmp, k[6], user_data);
	
	// Error
	OdeReal maxdiff = -std::numeric_limits<OdeReal>::infinity();
	for (size_t i = 0; i < N; i++) {
		OdeReal error = cur_dt * (
			+ (OdeReal)0.00123263888888888888888888888889 * k[0][i]
			- (OdeReal)0.00425277029050613956274333632824 * k[2][i]
			+ (OdeReal)0.03697916666666666666666666666667 * k[3][i]
			- (OdeReal)0.05086379716981132075471698113208 * k[4][i]
			+ (OdeReal)0.04190476190476190476190476190476 * k[5][i]
			- (OdeReal)0.025							  * k[6][i]);
		error = fabs(error);
		
		//OdeReal D = abstol[i] + reltol * (ytmp[i] + k[6][i] * cur_dt);
		OdeReal D = abstol + reltol * fabs(ytmp[i] + k[6][i] * cur_dt);
		OdeReal diff = error / D;
		maxdiff = (std::max)(maxdiff, diff);
	}

	return maxdiff;
}

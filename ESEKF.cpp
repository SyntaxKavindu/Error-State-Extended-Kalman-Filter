/*
 * ESEKF.cpp
 *
 *  Created on: Aug 16, 2026
 *      Author: KAVINDU
 */

#include "ESEKF.hpp"

#include <limits>
#include <type_traits>

// ---------------------------------------------------------------------------
// Build guards. These are not style preferences -- each one silently breaks a
// safety property of the filter, and each failure mode is invisible at runtime.
// ---------------------------------------------------------------------------
#if defined(__FAST_MATH__)
#error "ESEKF must not be built with -ffast-math / -ffinite-math-only: every divergence guard in this filter is a std::isfinite() check, and those flags license the compiler to assume NaN and Inf cannot occur and delete the guards. A diverged filter would then feed NaN straight to the control loop with hasDiverged() still returning false."
#endif

static_assert(sizeof(float) == 4, "ESEKF assumes IEEE-754 binary32 float");
static_assert(std::numeric_limits<float>::is_iec559,
              "ESEKF's guards rely on IEEE-754 semantics (NaN compares false, isfinite works)");
static_assert(std::is_trivially_destructible<ESEKF>::value,
              "ESEKF must be trivially destructible so a static/global instance needs no atexit slot");
static_assert(std::is_trivially_copyable<ESEKF>::value,
              "ESEKF must stay trivially copyable so a snapshot can be memcpy'd out for logging");

namespace
{
constexpr float kGravityMps2   = 9.80665f;
constexpr float kDefaultAttUnc = 0.05f;      // rad^2, initial attitude uncertainty
constexpr float kDefaultVelUnc = 0.25f;      // (m/s)^2
constexpr float kDefaultPosUnc = 1.0f;       // m^2
constexpr float kDefaultGyroBiasUnc  = 1.0e-4f; // (rad/s)^2
constexpr float kDefaultAccelBiasUnc = 1.0e-3f; // (m/s^2)^2
// WGS84 ellipsoid. The previous code used the equatorial radius for BOTH the
// north and the east conversion; the north conversion actually needs the
// meridional radius of curvature M, which is 0.67% smaller at the equator, so
// that shortcut put a ~0.7% scale error on the local frame's north axis.
constexpr double kWgs84A  = 6378137.0;          // semi-major axis, m
constexpr double kWgs84E2 = 6.69437999014e-3;   // first eccentricity squared

// Loose chi-square gate. dim(z) = 3 at 99.9% is 16.3, dim(z) = 1 at 99.9% is
// 10.8; 25 is deliberately permissive so the gate only ever catches a sample
// that is grossly inconsistent, never healthy transients.
constexpr float kDefaultNisGate = 25.0f;

// Sentinel for "this source has never been fused / rejected". Far enough in the
// past that every elapsed-time comparison reports a timeout. Double, to match
// the filter clock -- see the note on filter_time_ in the header.
constexpr double kNeverFused = -1000.0;

// Floor for a variance seeded by a state reset. See resetStateBlockCovariance().
constexpr float kMinSeedVariance = 1.0e-9f;

// How hard updateAccelerometer() de-weights a sample sitting at the motion
// gate: R is scaled by 1 + gain*(motion/gate)^2, so at the gate itself the
// measurement is charged sqrt(1 + 9) ~ 3.2x its nominal sigma. Chosen so a
// sample on the boundary still contributes something rather than falling off a
// cliff, while contributing far less than a clean at-rest reading.
constexpr float kAccelInflationGain = 9.0f;

// M_PI is a POSIX extension, not ISO C++; some arm-none-eabi/newlib setups hide
// it under a strict -std=c++20. Spell it out rather than depend on that.
constexpr double kPi     = 3.14159265358979323846;
constexpr double kTwoPi  = 6.28318530717958647692;
constexpr double kHalfPi = 1.57079632679489662;

// ---------------------------------------------------------------------------
// Input validation helpers.
//
// Every public setter runs its argument through these. That is not defensive
// clutter: on a real vehicle these values come from a parameter store, a
// ground-station link, or an auto-tune routine, and a single NaN written into
// R, Q or the gravity vector propagates into P on the next update and latches
// the filter into hasDiverged() -- in flight, with no way back except a reset
// and re-alignment. A rejected setter leaves the previous, known-good tuning
// in place, which is always the safer failure.
// ---------------------------------------------------------------------------
inline bool isFiniteVec(const Vector3f &v)
{
	return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}

inline bool isFiniteMat3(const float M[3][3])
{
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			if (!std::isfinite(M[i][j]))
			{
				return false;
			}
		}
	}
	return true;
}

// A measurement noise covariance must be finite and have a strictly positive
// diagonal. A zero or negative variance makes S = H P H^T + R singular or
// indefinite, and the update is then either dropped forever (matrixInverse3x3
// refuses) or fused with a gain that is not a Kalman gain at all.
inline bool isValidNoiseMat3(const float M[3][3])
{
	return isFiniteMat3(M) && M[0][0] > 0.0f && M[1][1] > 0.0f && M[2][2] > 0.0f;
}

// A process-noise or variance parameter: finite and non-negative. Zero is
// allowed (a state the caller asserts is noise-free) but must be deliberate.
inline bool isValidVariance(float v)
{
	return std::isfinite(v) && v >= 0.0f;
}

// Plausible geodetic angles in RADIANS. Latitude outside +-pi/2 is not a
// latitude, and the overwhelmingly likely cause is degrees being passed to a
// radians API -- the one mistake this filter's GPS path cannot detect any
// other way, and one that otherwise shows up only as a ~57x scale error on
// every horizontal position. Longitude is left wide (a caller may use either
// [-pi, pi] or [0, 2pi); wrapPi() handles both) but still bounded.
inline bool isValidLatLonRad(double lat, double lon)
{
	return std::isfinite(lat) && std::isfinite(lon) &&
	       lat >= -kHalfPi && lat <= kHalfPi &&
	       lon >= -kTwoPi && lon <= kTwoPi;
}

// Meridional (north) and prime-vertical (east) radii of curvature at `lat`.
// r_east already carries the cos(lat) factor, i.e. east = dlon * r_east.
inline void earthRadii(double lat, double &r_north, double &r_east)
{
	const double s  = std::sin(lat);
	const double w2 = 1.0 - kWgs84E2 * s * s;
	const double w  = std::sqrt(w2);
	r_north = kWgs84A * (1.0 - kWgs84E2) / (w2 * w);  // M
	r_east  = kWgs84A / w * std::cos(lat);            // N * cos(lat)
}

// Wraps a longitude difference into [-pi, pi] so a track crossing the
// antimeridian does not produce a 2*pi jump in the east coordinate. Branchless
// and bounded -- a while-loop here would hang on a garbage input.
inline double wrapPi(double a)
{
	if (a >= -kPi && a <= kPi)
	{
		return a; // the overwhelmingly common case
	}
	return a - kTwoPi * std::floor((a + kPi) / kTwoPi);
}

// Single-precision wrap, for the heading residual. Deliberately not a call
// through the double version: the Cortex-M4F/M7F FPU this targets is
// single-precision only, so every double operation is a software routine --
// tens of cycles each, on a path that runs at the magnetometer rate.
constexpr float kPiF    = 3.14159265358979f;
constexpr float kTwoPiF = 6.28318530717959f;

inline float wrapPiF(float a)
{
	if (a >= -kPiF && a <= kPiF)
	{
		return a; // the overwhelmingly common case: successive headings are close
	}
	return a - kTwoPiF * std::floor((a + kPiF) / kTwoPiF);
}

// Floor on the heading measurement variance, rad^2. sqrt(1e-6) is 1 mrad, or
// 0.06 deg: past that the magnetometer would be claiming a heading better than
// any hobby-grade part can deliver, and an over-confident yaw measurement pulls
// the gyro-bias states with it. Reached only if the caller sets an
// unrealistically small mag noise or the horizontal field reads very large.
constexpr float kMinYawVariance = 1.0e-6f;
}

// ---------------------------------------------------------------------------
// Constructor: bring every member to a well-defined, benign default so the
// filter is never run against uninitialized memory even if predict()/update()
// are called before inilialize().
// ---------------------------------------------------------------------------
ESEKF::ESEKF()
{
	q = {1.0f, 0.0f, 0.0f, 0.0f}; // identity orientation

	velocity_NED_ = {0.0f, 0.0f, 0.0f};
	position_NED_ = {0.0f, 0.0f, 0.0f};
	gyro_bias     = {0.0f, 0.0f, 0.0f};
	accel_bias    = {0.0f, 0.0f, 0.0f};

	std::memset(Q, 0, sizeof(Q));

	std::memset(R_accel, 0, sizeof(R_accel));
	std::memset(R_mag, 0, sizeof(R_mag));
	std::memset(R_gps, 0, sizeof(R_gps));
	std::memset(R_gps_vel, 0, sizeof(R_gps_vel));
	R_baro = 1.0f;
	accel_gate_threshold_ = 1.0f; // m/s^2 (~10% of g) -- tune with setAccelGateThreshold()
	accel_motion_         = 0.0f;
	accel_inflation_      = 1.0f;

	last_good_q_          = {1.0f, 0.0f, 0.0f, 0.0f};
	last_good_velocity_   = {0.0f, 0.0f, 0.0f};
	last_good_position_   = {0.0f, 0.0f, 0.0f};
	last_good_gyro_bias_  = {0.0f, 0.0f, 0.0f};
	last_good_accel_bias_ = {0.0f, 0.0f, 0.0f};
	have_last_good_       = false;
	fault_count_          = 0;

	// Generic MEMS-grade placeholders -- override with your IMU's actual
	// datasheet/Allan-variance figures via setImuNoiseParameters().
	gyro_noise_density_     = 0.005f;  // (rad/s)/sqrt(Hz)
	accel_noise_density_    = 0.02f;   // (m/s^2)/sqrt(Hz)
	gyro_bias_random_walk_  = 1.0e-4f; // (rad/s)/sqrt(s)
	accel_bias_random_walk_ = 1.0e-3f; // (m/s^2)/sqrt(s)

	g       = {0.0f, 0.0f, kGravityMps2}; // NED: down is positive z
	angular_rate_ = {0.0f, 0.0f, 0.0f};
	mag_ref_ = {1.0f, 0.0f, 0.0f};        // placeholder until inilialize() sets it
	gps_ref_ = {0.0f, 0.0f, 0.0f};        // placeholder until inilialize() sets it

	gps_ref_lat_rad_ = 0.0;
	gps_ref_lon_rad_ = 0.0;

	initialized_ = false;
	diverged_    = false;
	dt_last_     = 0.0f;
	altitude_    = 0.0f;
	baro_ref_    = 0.0f;
	nis_gate_    = kDefaultNisGate;
	mag_declination_ = 0.0f;
	glitch_radius_   = 25.0f; // EKF3 EK3_GLITCH_RAD default

	last_pos_reset_delta_ = {0.0f, 0.0f, 0.0f};
	last_vel_reset_delta_ = {0.0f, 0.0f, 0.0f};
	last_pos_reset_t_     = kNeverFused;
	last_vel_reset_t_     = kNeverFused;
	last_pos_reset_exceeded_glitch_ = false;

	filter_time_           = 0.0;
	last_gps_pos_fuse_t_   = kNeverFused;
	last_gps_vel_fuse_t_   = kNeverFused;
	last_baro_fuse_t_      = kNeverFused;
	last_mag_fuse_t_       = kNeverFused;
	last_gps_pos_reject_t_ = kNeverFused;
	last_gps_vel_reject_t_ = kNeverFused;
	last_baro_reject_t_    = kNeverFused;
	last_mag_reject_t_     = kNeverFused;

	pos_test_ratio_ = 0.0f;
	vel_test_ratio_ = 0.0f;
	hgt_test_ratio_ = 0.0f;
	mag_test_ratio_ = 0.0f;

	last_test_ratio_   = 0.0f;
	last_update_gated_ = false;

	// Reasonable default measurement noise / position process noise so the
	// filter is usable even if the caller never calls set*Noise(). Unlike P,
	// these are tuning parameters, not filter state, so inilialize() does not
	// reset them on a re-init.
	for (int i = 0; i < 3; ++i)
	{
		Q[6 + i][6 + i] = 1.0e-6f; // small regularizing term on the position block

		R_accel[i][i] = 0.05f;
		R_mag[i][i]   = 0.05f;
		R_gps[i][i]   = 1.0f;
	}

	// GPS velocity defaults follow EKF3's EK3_VELNE_M_NSE / EK3_VELD_M_NSE
	// (0.3 m/s horizontal, 0.5 m/s vertical for copter).
	setGPSVelocityNoiseSigma(0.3f, 0.5f);

	// Default IMU noise-density figures (typical low-cost MEMS IMU, e.g.
	// MPU6050-class). Overwrite with your actual datasheet values via
	// setImuNoiseParameters() -- this only fills the attitude/velocity/bias
	// blocks of Q as a reasonable starting point.
	setImuNoiseParameters(
		1.0e-3f,  // gyro_noise_density,  (rad/s)/sqrt(Hz)
		2.0e-2f,  // accel_noise_density, (m/s^2)/sqrt(Hz)
		3.0e-4f,  // gyro_bias_random_walk,  (rad/s)/sqrt(s)
		3.0e-3f); // accel_bias_random_walk, (m/s^2)/sqrt(s)

	resetCovarianceDefaults();
}

// ---------------------------------------------------------------------------
// resetCovarianceDefaults(): sets P to its default diagonal. Shared by the
// constructor and inilialize() so the two never drift out of sync.
// ---------------------------------------------------------------------------
void ESEKF::resetCovarianceDefaults()
{
	std::memset(P, 0, sizeof(P));

	for (int i = 0; i < 3; ++i)
	{
		P[i][i]           = kDefaultAttUnc;
		P[3 + i][3 + i]   = kDefaultVelUnc;
		P[6 + i][6 + i]   = kDefaultPosUnc;
		P[9 + i][9 + i]   = kDefaultGyroBiasUnc;
		P[12 + i][12 + i] = kDefaultAccelBiasUnc;
	}
}

// ---------------------------------------------------------------------------
// Small math helpers that inilialize() depends on.
// ---------------------------------------------------------------------------
float ESEKF::vectorNorm(const Vector3f &a) const
{
	return std::sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
}

Vector3f ESEKF::crossProduct(const Vector3f &a, const Vector3f &b) const
{
	return Vector3f{
		a.y * b.z - a.z * b.y,
		a.z * b.x - a.x * b.z,
		a.x * b.y - a.y * b.x
	};
}

Quaternionf ESEKF::quaternionNormalize(const Quaternionf &q_in) const
{
	float n = std::sqrt(q_in.w * q_in.w + q_in.x * q_in.x +
	                     q_in.y * q_in.y + q_in.z * q_in.z);

	if (n < 1.0e-9f)
	{
		return Quaternionf{1.0f, 0.0f, 0.0f, 0.0f};
	}

	float inv_n = 1.0f / n;
	return Quaternionf{q_in.w * inv_n, q_in.x * inv_n, q_in.y * inv_n, q_in.z * inv_n};
}

Quaternionf ESEKF::quaternionMultiply(const Quaternionf &q1, const Quaternionf &q2) const
{
	return Quaternionf{
		q1.w * q2.w - q1.x * q2.x - q1.y * q2.y - q1.z * q2.z,
		q1.w * q2.x + q1.x * q2.w + q1.y * q2.z - q1.z * q2.y,
		q1.w * q2.y - q1.x * q2.z + q1.y * q2.w + q1.z * q2.x,
		q1.w * q2.z + q1.x * q2.y - q1.y * q2.x + q1.z * q2.w
	};
}

Quaternionf ESEKF::quaternionConjugate(const Quaternionf &q_in) const
{
	return Quaternionf{q_in.w, -q_in.x, -q_in.y, -q_in.z};
}

// Roll/pitch/yaw (rad, body 3-2-1 Euler) -> orientation quaternion (body -> NED)
Quaternionf ESEKF::quaternionFromEuler(float roll, float pitch, float yaw) const
{
	float cr = std::cos(roll * 0.5f),  sr = std::sin(roll * 0.5f);
	float cp = std::cos(pitch * 0.5f), sp = std::sin(pitch * 0.5f);
	float cy = std::cos(yaw * 0.5f),   sy = std::sin(yaw * 0.5f);

	Quaternionf out;
	out.w = cr * cp * cy + sr * sp * sy;
	out.x = sr * cp * cy - cr * sp * sy;
	out.y = cr * sp * cy + sr * cp * sy;
	out.z = cr * cp * sy - sr * sp * cy;

	return quaternionNormalize(out);
}

// Rotates v from body frame to NED frame using q (body -> NED).
Vector3f ESEKF::rotateVector(const Quaternionf &q_in, const Vector3f &v) const
{
	Quaternionf v_quat{0.0f, v.x, v.y, v.z};
	Quaternionf q_conj = quaternionConjugate(q_in);
	Quaternionf result = quaternionMultiply(quaternionMultiply(q_in, v_quat), q_conj);

	return Vector3f{result.x, result.y, result.z};
}

// ---------------------------------------------------------------------------
// inilialize(): compute an initial attitude from the accelerometer (leveling)
// and magnetometer (tilt-compensated heading), seed position from GPS and
// altitude, zero the velocity, and reset the covariance to its defaults.
// ---------------------------------------------------------------------------
void ESEKF::initialize(const Vector3f &accel, const Vector3f &mag, float altitude, const Vector3f &gps)
{
	// --- Reject an unusable alignment BEFORE touching any state ---
	//
	// Alignment is the one moment the filter takes a measurement entirely on
	// trust: there is no prior to gate it against, so whatever comes in here
	// becomes the attitude, the origin and the datum. A bad sample at this
	// point does not produce a bad estimate that converges away, it produces a
	// filter that is confidently wrong for the whole flight.
	//
	// Failing leaves isInitialized() false, and every predict()/update() is
	// then a no-op until the caller retries with good data. CHECK IT.
	initialized_ = false;

	// Drop the rollback point: it belongs to the previous alignment, and a state
	// from before a re-alignment is not something this run may fall back on. It
	// is re-armed by the health check at the bottom of this function, which is
	// also what makes a garbage alignment latch diverged_ rather than repair
	// itself out of a state that no longer applies.
	have_last_good_ = false;

	if (!isFiniteVec(accel) || !isFiniteVec(mag) || !std::isfinite(altitude) ||
	    !std::isfinite(gps.z) || !isValidLatLonRad((double)gps.x, (double)gps.y))
	{
		return; // non-finite sample, or lat/lon that cannot be radians (degrees?)
	}

	// --- Leveling from accelerometer (assumes the vehicle is stationary) ---
	//
	// A stationary FRD accelerometer measures specific force R^T * (-g_ned) with
	// g_ned = [0,0,+9.81], i.e. a = [ g*sin(pitch), -g*sin(roll)cos(pitch),
	// -g*cos(roll)cos(pitch) ]. Inverting that gives the two atan2 forms below;
	// both are exact for NED + FRD (verified against a full roll/pitch/yaw grid).
	const float accel_norm = vectorNorm(accel);

	// The reading must actually look like gravity.
	//
	// Alignment takes this sample entirely on trust -- there is no prior to gate
	// it against -- so a dead, saturated or stuck accelerometer aligns the filter
	// just as confidently as a good one. A dead part reading zero hands back
	// roll = pitch = 0, which is indistinguishable from a vehicle that happens to
	// be level; a saturated one reading a huge constant gives a perfectly
	// plausible attitude from a direction that means nothing. Both are then
	// unrecoverable, because the accelerometer update's own motion gate rejects
	// exactly the samples that would have corrected the error.
	//
	// The band is deliberately loose -- half of local gravity either way -- so it
	// never refuses a real stationary reading from an uncalibrated part or a
	// vehicle being held by hand. It is here to catch a broken sensor, not to
	// judge stillness: a caller that wants a tight stillness gate should apply
	// one before calling, which is what the vehicle layer does.
	//
	// Compared against |g| rather than a literal 9.81 so it stays correct for a
	// caller working in g-units or on another body, via setGravity().
	const float gravity_mag = vectorNorm(g);
	if (!(gravity_mag > 1.0e-6f) ||
	    std::fabs(accel_norm - gravity_mag) > 0.5f * gravity_mag)
	{
		return;
	}

	const float roll  = std::atan2(-accel.y, -accel.z);
	const float pitch = std::atan2(accel.x, std::sqrt(accel.y * accel.y + accel.z * accel.z));

	// --- Tilt-compensated magnetic heading ---
	//
	// With R_b^n = Rz(yaw) * Ry(pitch) * Rx(roll), levelling the body-frame field
	// requires m_level = Ry(pitch) * Rx(roll) * m_body -- in THAT order. The
	// widely copied "Xh = mx*cp + mz*sp / Yh = mx*sr*sp + my*cr - mz*sr*cp" pair
	// is Rx(roll) * Ry(pitch) * m_body instead, i.e. the two rotations applied in
	// the wrong order. It happens to be exact whenever roll or pitch is zero,
	// which is why it survives casual testing, but it is wrong whenever both are
	// non-zero: 3.6 deg of heading error at 20 deg/20 deg, and up to 34 deg at
	// large tilt. The correct levelling is used here.
	float cr = std::cos(roll),  sr = std::sin(roll);
	float cp = std::cos(pitch), sp = std::sin(pitch);

	// Ry(pitch) * Rx(roll) applied to mag:
	//   row 0 = [ cp,  sp*sr,  sp*cr ]
	//   row 1 = [  0,     cr,    -sr ]
	float mx = mag.x * cp + mag.y * sp * sr + mag.z * sp * cr;
	float my =              mag.y * cr      - mag.z * sr;

	// atan2(-my, mx) is the MAGNETIC heading. Adding the declination (angle from
	// true north to magnetic north, positive east) converts it to TRUE heading.
	// EKF3 forms the same quantity as
	//   yawAngMeasured = wrap_PI(-atan2(magNED.y, magNED.x) + MagDeclination()).
	// A zero horizontal field carries no heading. That is not necessarily an
	// error -- a vehicle with no magnetometer may legitimately pass {0,0,0} --
	// so yaw starts at zero and the local frame is simply not north-aligned
	// until some other source (GPS velocity, an external heading) constrains
	// it. It is NOT a reason to refuse alignment, unlike a dead accelerometer:
	// roll and pitch remain fully observable without a magnetometer.
	float yaw = 0.0f;
	if (!(mx == 0.0f && my == 0.0f))
	{
		yaw = std::atan2(-my, mx) + mag_declination_;
	}

	q = quaternionFromEuler(roll, pitch, yaw);

	// --- Position / velocity ---
	// gps is a raw reading (lat/lon/alt), not a local NED position, so it cannot
	// be assigned to position_NED_ directly. Store it as the origin of the local
	// NED frame; the filter starts at that origin by definition.
	setGPSReferencePrecise((double)gps.x, (double)gps.y, gps.z);
	position_NED_   = {0.0f, 0.0f, 0.0f}; // local NED origin -- absolute altitude is reconstructed via gps_ref_.z (see getAltitude())
	altitude_       = altitude;

	// --- Barometric datum ---
	// The barometer's own first reading is the reference the barometer residual
	// is formed against. It is NOT gps_ref_.z: barometric altitude (whatever
	// pressure datum the driver uses) and GPS altitude (ellipsoidal or geoid
	// height) differ by a large, weather- and datum-dependent offset. Mixing the
	// two puts that whole offset straight into the vertical position residual.
	baro_ref_ = altitude;

	velocity_NED_ = {0.0f, 0.0f, 0.0f};

	// --- Biases start unknown ---
	gyro_bias  = {0.0f, 0.0f, 0.0f};
	accel_bias = {0.0f, 0.0f, 0.0f};

	angular_rate_ = {0.0f, 0.0f, 0.0f};

	// A reset reported against the PREVIOUS alignment says nothing about this
	// one, and the clock it was stamped on is about to restart at zero -- which
	// would make an old reset look like it happened at t = 0 of the new run and
	// raise position_reset for the first seconds of flight.
	last_pos_reset_delta_ = {0.0f, 0.0f, 0.0f};
	last_vel_reset_delta_ = {0.0f, 0.0f, 0.0f};
	last_pos_reset_t_     = kNeverFused;
	last_vel_reset_t_     = kNeverFused;
	last_pos_reset_exceeded_glitch_ = false;

	// --- Magnetic reference: this is what makes yaw absolute or relative ---
	//
	// The obvious construction, mag_ref_ = rotateVector(q, mag), rotates the
	// first sample through the attitude that was just derived FROM that same
	// sample. The result is self-consistent by construction: the very first
	// magnetometer innovation is identically zero, and it stays zero for any
	// initial yaw whatsoever. The magnetometer then only holds yaw wherever
	// initialization happened to put it and supplies no absolute heading -- an
	// initial heading error becomes permanent and unobservable.
	//
	// Instead, build the NED reference field so its HORIZONTAL component points
	// at true magnetic north, using the declination. Only the levelled (roll and
	// pitch) rotation is used, which is trustworthy because it comes from
	// gravity rather than from the magnetometer; the horizontal magnitude and
	// the down component are taken from the measurement so the reference keeps
	// the sensor's own scale and local inclination. This mirrors EKF3's
	// alignMagStateDeclination(), which rotates the NE earth-field states so the
	// declination matches the published World Magnetic Model value.
	//
	// (mx, my) computed above IS the levelled field, so its magnitude is the
	// horizontal field strength. The down component follows from the same
	// levelling: row 2 of Ry(pitch)*Rx(roll).
	{
		const float mh = std::sqrt(mx * mx + my * my);
		const float md = -mag.x * sp + mag.y * cp * sr + mag.z * cp * cr;

		mag_ref_ = Vector3f{mh * std::cos(mag_declination_),
		                    mh * std::sin(mag_declination_),
		                    md};
	}

	// --- Reset covariance to initial uncertainty ---
	resetCovarianceDefaults();

	dt_last_     = 0.0f;
	diverged_    = false;

	// Restart the aiding clock: nothing has been fused against this alignment.
	filter_time_           = 0.0;
	last_gps_pos_fuse_t_   = kNeverFused;
	last_gps_vel_fuse_t_   = kNeverFused;
	last_baro_fuse_t_      = kNeverFused;
	last_mag_fuse_t_       = kNeverFused;
	last_gps_pos_reject_t_ = kNeverFused;
	last_gps_vel_reject_t_ = kNeverFused;
	last_baro_reject_t_    = kNeverFused;
	last_mag_reject_t_     = kNeverFused;
	pos_test_ratio_ = 0.0f;
	vel_test_ratio_ = 0.0f;
	hgt_test_ratio_ = 0.0f;
	mag_test_ratio_ = 0.0f;

	initialized_ = true;

	// If the caller handed us garbage, fail loudly here rather than silently
	// running a NaN filter for the rest of the flight.
	if (!checkFinite())
	{
		initialized_ = false;
	}
}

// ---------------------------------------------------------------------------
// Small vector helpers.
// ---------------------------------------------------------------------------
Vector3f ESEKF::vectorAdd(const Vector3f &a, const Vector3f &b) const
{
	return Vector3f{a.x + b.x, a.y + b.y, a.z + b.z};
}

Vector3f ESEKF::vectorSub(const Vector3f &a, const Vector3f &b) const
{
	return Vector3f{a.x - b.x, a.y - b.y, a.z - b.z};
}

Vector3f ESEKF::vectorScale(const Vector3f &a, float s) const
{
	return Vector3f{a.x * s, a.y * s, a.z * s};
}

// 3x3 skew-symmetric ("cross-product") matrix of v, so that skew(v) * x == v cross x.
void ESEKF::skewSymmetric(const Vector3f &v, float M_out[3][3]) const
{
	M_out[0][0] =  0.0f;   M_out[0][1] = -v.z;   M_out[0][2] =  v.y;
	M_out[1][0] =  v.z;    M_out[1][1] =  0.0f;  M_out[1][2] = -v.x;
	M_out[2][0] = -v.y;    M_out[2][1] =  v.x;   M_out[2][2] =  0.0f;
}

// Exact exponential map: rotation vector delta_theta -> unit quaternion.
//   dq = [cos(|dt|/2), sin(|dt|/2) * dt/|dt|]
//
// This replaces the previous first-order form dq = normalize([1, dtheta/2]).
// That form represents a rotation of 2*atan(|dtheta|/2) rather than |dtheta|,
// i.e. it under-rotates by |dtheta|^3/12 every step. Integrated at rate omega
// and step dt that is a systematic attitude drift of omega^3 * dt^2 / 12 rad/s
// -- 0.12 deg/s at 10 rad/s and 200 Hz, 0.95 deg/s at 20 rad/s -- which is a
// RATE-DEPENDENT error that the gyro-bias state cannot absorb. It also matters
// for error-state injection, where a large correction (e.g. after a GPS outage)
// is no longer a small angle at all.
//
// Cost on a Cortex-M4F is one sqrtf plus one sincos, a few hundred nanoseconds.
Quaternionf ESEKF::quaternionFromRotationVector(const Vector3f &delta_theta) const
{
	const float a2 = delta_theta.x * delta_theta.x +
	                 delta_theta.y * delta_theta.y +
	                 delta_theta.z * delta_theta.z;
	const float a  = std::sqrt(a2);

	if (a < 1.0e-8f)
	{
		// Below this the series and the closed form agree to well within float
		// precision, and the closed form would divide by ~0.
		return quaternionNormalize(
			Quaternionf{1.0f, 0.5f * delta_theta.x, 0.5f * delta_theta.y, 0.5f * delta_theta.z});
	}

	const float half = 0.5f * a;
	const float sc   = std::sin(half) / a;

	return Quaternionf{std::cos(half), sc * delta_theta.x, sc * delta_theta.y, sc * delta_theta.z};
}

// Body -> NED rotation matrix from the orientation quaternion.
void ESEKF::quaternionToRotationMatrix(const Quaternionf &q_in, float R_out[3][3]) const
{
	float w = q_in.w, x = q_in.x, y = q_in.y, z = q_in.z;

	R_out[0][0] = 1 - 2 * (y * y + z * z);
	R_out[0][1] = 2 * (x * y - w * z);
	R_out[0][2] = 2 * (x * z + w * y);

	R_out[1][0] = 2 * (x * y + w * z);
	R_out[1][1] = 1 - 2 * (x * x + z * z);
	R_out[1][2] = 2 * (y * z - w * x);

	R_out[2][0] = 2 * (x * z - w * y);
	R_out[2][1] = 2 * (y * z + w * x);
	R_out[2][2] = 1 - 2 * (x * x + y * y);
}

// ---------------------------------------------------------------------------
// The only fixed-size matrix product the filter needs. Fully unrolled inner
// loop -- at 3x3 the loop overhead is comparable to the arithmetic.
// ---------------------------------------------------------------------------
void ESEKF::matrixMultiply3x3(const float A[3][3], const float B[3][3], float C_out[3][3]) const
{
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			C_out[i][j] = A[i][0] * B[0][j] + A[i][1] * B[1][j] + A[i][2] * B[2][j];
		}
	}
}

// ---------------------------------------------------------------------------
// predictCovariance(): P <- F P F^T + Q*dt, WITHOUT ever materializing F.
//
// F = I + A*dt is the discrete error-state transition. Written as blocks over
// [dtheta, dv, dp, dbg, dba], with everything evaluated at the mid-interval
// attitude R (the same linearization point predict() integrates velocity at):
//
//        | I - [w]x dt      0      0   -I dt      0   |
//        |  -R[f]x dt       I      0      0    -R dt  |
//   F =  |      0         I dt     I      0       0   |
//        |      0           0      0      I       0   |
//        |      0           0      0      0       I   |
//
// 195 of F's 225 entries are structurally zero and 15 more are the identity,
// so a dense F P F^T spends 6750 multiply-accumulates computing with zeros and
// needs 1.8 kB of stack (900 bytes for F, 900 for the F*P temporary). That is
// not a micro-optimization on an MCU: an STM32 FreeRTOS task is routinely
// given a 1-2 kB stack, and predict() is the highest-rate entry point in the
// filter, so the dense form is a stack overflow in the hot path.
//
// Left-multiplying by F is a set of row operations, right-multiplying by F^T
// the mirror-image column operations, and each one touches only the three rows
// (columns) of one state block:
//
//   rows 0..2 :  theta <- theta - dt*[w]x theta - dt*bg
//   rows 3..5 :  v     <- v     - dt*R[f]x theta - dt*R ba
//   rows 6..8 :  p     <- p     + dt*v
//   rows 9..14:  unchanged (both biases are random walks)
//
// ORDER MATTERS and is the one subtle thing here. Each block must still see
// the ORIGINAL value of every block it reads, so they are applied 6..8 first
// (it consumes rows 3..5), then 3..5 (consumes rows 0..2 and 12..14), then
// 0..2 (consumes itself -- hence the t0/t1/t2 snapshot -- and rows 9..11).
// Applying them in the natural 0..14 order would feed already-updated rows
// into later blocks and silently compute F P with the wrong F.
//
// Result: exactly the same P as the dense product (verified numerically
// against an independent dense implementation to float round-off), for ~1000
// multiply-accumulates and a few dozen bytes of stack.
//
// Note on Q: the noise-input matrix G is the identity here, which is correct
// only because every Q block is ISOTROPIC (sigma^2 * I). The accelerometer
// noise actually enters the error dynamics as -R*n_a, so its contribution is
// R*(sigma_a^2 I)*R^T = sigma_a^2 I -- rotation-invariant. If you ever make a Q
// block anisotropic you must reinstate G and rotate it.
// ---------------------------------------------------------------------------
void ESEKF::predictCovariance(const Vector3f &omega, const Vector3f &f,
                               const float R[3][3], float dt)
{
	float S_omega[3][3];
	skewSymmetric(omega, S_omega);

	float S_f[3][3];
	skewSymmetric(f, S_f);

	float RSf[3][3];
	matrixMultiply3x3(R, S_f, RSf);

	// ---- Pass 1: P <- F P, as row operations, one column at a time. ----
	// Row operations combine rows within a single column, so the columns are
	// independent and the whole pass is safe to do in place.
	for (int j = 0; j < ESEKF_STATE_DIM; ++j)
	{
		// Snapshot the attitude rows: rows 0..2 are read by the velocity block
		// below AND by their own update, which overwrites them.
		const float t0 = P[0][j], t1 = P[1][j], t2 = P[2][j];
		const float bg[3] = {P[9][j],  P[10][j], P[11][j]};
		const float ba[3] = {P[12][j], P[13][j], P[14][j]};

		// p <- p + dt*v   (must precede the v update)
		for (int i = 0; i < 3; ++i)
		{
			P[6 + i][j] += dt * P[3 + i][j];
		}

		// v <- v - dt*R[f]x theta - dt*R ba
		for (int i = 0; i < 3; ++i)
		{
			P[3 + i][j] -= dt * (RSf[i][0] * t0 + RSf[i][1] * t1 + RSf[i][2] * t2)
			             + dt * (R[i][0] * ba[0] + R[i][1] * ba[1] + R[i][2] * ba[2]);
		}

		// theta <- theta - dt*[w]x theta - dt*bg
		for (int i = 0; i < 3; ++i)
		{
			P[i][j] -= dt * (S_omega[i][0] * t0 + S_omega[i][1] * t1 + S_omega[i][2] * t2)
			         + dt * bg[i];
		}
	}

	// ---- Pass 2: P <- P F^T, the mirror image, one row at a time. ----
	for (int i = 0; i < ESEKF_STATE_DIM; ++i)
	{
		const float t0 = P[i][0], t1 = P[i][1], t2 = P[i][2];
		const float bg[3] = {P[i][9],  P[i][10], P[i][11]};
		const float ba[3] = {P[i][12], P[i][13], P[i][14]};

		for (int j = 0; j < 3; ++j)
		{
			P[i][6 + j] += dt * P[i][3 + j];
		}

		for (int j = 0; j < 3; ++j)
		{
			P[i][3 + j] -= dt * (RSf[j][0] * t0 + RSf[j][1] * t1 + RSf[j][2] * t2)
			             + dt * (R[j][0] * ba[0] + R[j][1] * ba[1] + R[j][2] * ba[2]);
		}

		for (int j = 0; j < 3; ++j)
		{
			P[i][j] -= dt * (S_omega[j][0] * t0 + S_omega[j][1] * t1 + S_omega[j][2] * t2)
			         + dt * bg[j];
		}
	}

	// ---- Process noise ----
	for (int i = 0; i < ESEKF_STATE_DIM; ++i)
	{
		for (int j = 0; j < ESEKF_STATE_DIM; ++j)
		{
			P[i][j] += Q[i][j] * dt;
		}
	}

	symmetrizeCovariance();
}

// ---------------------------------------------------------------------------
// predict(): propagates the nominal state (strapdown mechanization) and the
// error-state covariance forward by dt using gyro/accel measurements.
// ---------------------------------------------------------------------------
void ESEKF::predict(const Vector3f &gyro, const Vector3f &accel, float dt)
{
	if (!initialized_ || diverged_)
	{
		return;
	}

	// `dt <= 0.0f` alone does NOT reject a NaN: every comparison against NaN is
	// false in IEEE-754, so a NaN dt walks straight through such a guard and
	// turns the entire state and covariance into NaN with no recovery path.
	// Test for finiteness explicitly, and clamp a too-large step rather than
	// letting a scheduling stall integrate a huge dt.
	if (!std::isfinite(dt) || dt <= 0.0f)
	{
		return;
	}
	if (!std::isfinite(gyro.x)  || !std::isfinite(gyro.y)  || !std::isfinite(gyro.z) ||
	    !std::isfinite(accel.x) || !std::isfinite(accel.y) || !std::isfinite(accel.z))
	{
		return;
	}
	if (dt > ESEKF_MAX_PREDICT_DT)
	{
		dt = ESEKF_MAX_PREDICT_DT;
	}

	Vector3f omega = vectorSub(gyro, gyro_bias);  // bias-corrected angular rate
	Vector3f f     = vectorSub(accel, accel_bias); // bias-corrected specific force

	angular_rate_ = omega; // cache for getAngularRate() -- rate controllers read this

	// --- Attitude propagation ---
	// q_k = q_{k-1} (x) Exp(omega*dt). Right (body-frame) multiplication, which
	// is the same side the error state is defined on -- see injectErrorState().
	Vector3f delta_theta = vectorScale(omega, dt);
	Quaternionf q_prev   = q;
	q = quaternionNormalize(quaternionMultiply(q_prev, quaternionFromRotationVector(delta_theta)));

	// --- Mid-interval attitude ---
	// Specific force must be rotated by the attitude that held DURING the
	// interval, not at either endpoint. Rotating with q_{k-1} (forward Euler) or
	// with q_k (backward Euler) are both first-order accurate but each carries a
	// half-step of attitude misalignment, worth R*(omega x f)*dt^2/2 of velocity
	// error per step -- a systematic, coherent error under sustained rotation.
	// The midpoint attitude makes the velocity integration second-order accurate
	// for the price of one extra quaternion product, and it is also the correct
	// linearization point for F below.
	Quaternionf q_mid = quaternionNormalize(
		quaternionMultiply(q_prev, quaternionFromRotationVector(vectorScale(omega, 0.5f * dt))));

	float R_mid[3][3];
	quaternionToRotationMatrix(q_mid, R_mid);

	// --- Velocity propagation (specific force -> NED, plus gravity) ---
	// v_dot = R_b^n * f + g, with g = [0,0,+9.80665] in NED (down positive).
	// At rest f = R^T*(-g) so R*f = -g and v_dot = 0, as it must be.
	Vector3f f_ned{
		R_mid[0][0] * f.x + R_mid[0][1] * f.y + R_mid[0][2] * f.z,
		R_mid[1][0] * f.x + R_mid[1][1] * f.y + R_mid[1][2] * f.z,
		R_mid[2][0] * f.x + R_mid[2][1] * f.y + R_mid[2][2] * f.z};
	Vector3f accel_ned = vectorAdd(f_ned, g); // true acceleration in NED

	Vector3f velocity_prev = velocity_NED_;
	velocity_NED_ = vectorAdd(velocity_NED_, vectorScale(accel_ned, dt));

	// --- Position propagation (trapezoidal, using the average of old/new velocity) ---
	Vector3f velocity_avg = vectorScale(vectorAdd(velocity_prev, velocity_NED_), 0.5f);
	position_NED_ = vectorAdd(position_NED_, vectorScale(velocity_avg, dt));

	// gyro_bias / accel_bias: modeled as a random walk, so the nominal value
	// itself does not change in predict() -- only its uncertainty grows, via Q.

	// --- Error-state covariance propagation ---
	// F is applied in structured form, never built -- see predictCovariance().
	predictCovariance(omega, f, R_mid, dt);

	dt_last_ = dt;
	filter_time_ += (double)dt; // monotonic filter clock, drives every aiding timeout
	checkFinite();
}

// ---------------------------------------------------------------------------
// eulerFromQuaternion(): 3-2-1 (roll, pitch, yaw) Euler angles, returned as
// (x = roll, y = pitch, z = yaw), radians.
// ---------------------------------------------------------------------------
Vector3f ESEKF::eulerFromQuaternion(const Quaternionf &q_in) const
{
	Vector3f euler;

	// roll (x-axis rotation)
	float sinr_cosp = 2.0f * (q_in.w * q_in.x + q_in.y * q_in.z);
	float cosr_cosp = 1.0f - 2.0f * (q_in.x * q_in.x + q_in.y * q_in.y);
	euler.x = std::atan2(sinr_cosp, cosr_cosp);

	// pitch (y-axis rotation), clamp to avoid NaN from asin at the poles
	float sinp = 2.0f * (q_in.w * q_in.y - q_in.z * q_in.x);
	if (sinp > 1.0f)  sinp = 1.0f;
	if (sinp < -1.0f) sinp = -1.0f;
	euler.y = std::asin(sinp);

	// yaw (z-axis rotation)
	float siny_cosp = 2.0f * (q_in.w * q_in.z + q_in.x * q_in.y);
	float cosy_cosp = 1.0f - 2.0f * (q_in.y * q_in.y + q_in.z * q_in.z);
	euler.z = std::atan2(siny_cosp, cosy_cosp);

	return euler;
}

// ---------------------------------------------------------------------------
// matrixInverse3x3(): closed-form inverse via the adjugate / cofactor method.
// Returns false (leaving A_inv_out untouched) if A is singular.
//
// SINGULARITY TEST -- this must be RELATIVE, not absolute.
//
// The previous test was `|det| < 1e-12f`. A determinant of a 3x3 has the cube
// of the matrix's units, so that constant silently encodes an assumption about
// what the measurement is measured in. It is harmless for the accelerometer
// (|a| ~ 9.8 m/s^2, det(S) ~ 1e-2) and catastrophic for the magnetometer.
//
// For a magnetometer in Gauss, |mag| ~ 0.47 and R_mag ~ 1.6e-5, and
// H = skew(v) is rank 2 (v lies in its own null space), so one eigenvalue of
// S = H P H^T + R is essentially just R. det(S) ~ R * lambda1 * lambda2 lands
// around 1e-13 -- BELOW the cutoff -- as soon as the attitude covariance
// converges. Measured: the magnetometer update fuses normally while attitude
// sigma is above ~1.8 deg and is then silently dropped forever below it. S is
// never actually singular here (R is positive definite, so S is too); the
// guard was simply firing on a healthy matrix.
//
// The failure is invisible: the magnetometer is the only yaw reference, so yaw
// quietly degrades to unaided gyro dead-reckoning at exactly the point the
// filter appears to have converged, and every getter keeps returning plausible
// numbers.
//
// The test below compares |det| against the cube of the matrix's own scale, so
// it is dimensionless and rejects only genuine rank deficiency.
// ---------------------------------------------------------------------------
bool ESEKF::matrixInverse3x3(const float A[3][3], float A_inv_out[3][3]) const
{
	float det = A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1])
	          - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
	          + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);

	// Scale of the matrix: the largest magnitude entry. det of a well
	// conditioned matrix is O(scale^3), so |det| / scale^3 is dimensionless and
	// is roughly the product of the normalized singular values.
	float scale = 0.0f;
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			const float a = std::fabs(A[i][j]);
			if (a > scale)
			{
				scale = a;
			}
		}
	}

	// 1e-6 is ~10x float32's relative precision: below that the inverse would
	// be dominated by round-off. The absolute floor catches an all-zero matrix.
	const float scale3 = scale * scale * scale;
	if (!(scale > 0.0f) || !(std::fabs(det) > 1.0e-6f * scale3) || !std::isfinite(det))
	{
		return false; // genuinely singular / ill-conditioned, refuse to invert
	}

	float inv_det = 1.0f / det;

	A_inv_out[0][0] =  (A[1][1] * A[2][2] - A[1][2] * A[2][1]) * inv_det;
	A_inv_out[0][1] = -(A[0][1] * A[2][2] - A[0][2] * A[2][1]) * inv_det;
	A_inv_out[0][2] =  (A[0][1] * A[1][2] - A[0][2] * A[1][1]) * inv_det;

	A_inv_out[1][0] = -(A[1][0] * A[2][2] - A[1][2] * A[2][0]) * inv_det;
	A_inv_out[1][1] =  (A[0][0] * A[2][2] - A[0][2] * A[2][0]) * inv_det;
	A_inv_out[1][2] = -(A[0][0] * A[1][2] - A[0][2] * A[1][0]) * inv_det;

	A_inv_out[2][0] =  (A[1][0] * A[2][1] - A[1][1] * A[2][0]) * inv_det;
	A_inv_out[2][1] = -(A[0][0] * A[2][1] - A[0][1] * A[2][0]) * inv_det;
	A_inv_out[2][2] =  (A[0][0] * A[1][1] - A[0][1] * A[1][0]) * inv_det;

	return true;
}

// ---------------------------------------------------------------------------
// injectErrorState(): folds a 15x1 error state into the nominal state.
// dx layout: [0:3]=dtheta, [3:6]=dv, [6:9]=dp, [9:12]=dbg, [12:15]=dba.
// ---------------------------------------------------------------------------
void ESEKF::injectErrorState(const float dx[ESEKF_STATE_DIM])
{
	// Public entry point, so the correction is not necessarily one this filter
	// produced. A single non-finite element would be spread across the whole
	// nominal state by the quaternion product below.
	for (int i = 0; i < ESEKF_STATE_DIM; ++i)
	{
		if (!std::isfinite(dx[i]))
		{
			return;
		}
	}

	Vector3f dtheta{dx[0], dx[1], dx[2]};
	Quaternionf dq = quaternionFromRotationVector(dtheta);
	q = quaternionNormalize(quaternionMultiply(q, dq)); // RIGHT/local multiplication

	velocity_NED_ = vectorAdd(velocity_NED_, Vector3f{dx[3], dx[4], dx[5]});
	position_NED_ = vectorAdd(position_NED_, Vector3f{dx[6], dx[7], dx[8]});
	gyro_bias     = vectorAdd(gyro_bias, Vector3f{dx[9], dx[10], dx[11]});
	accel_bias    = vectorAdd(accel_bias, Vector3f{dx[12], dx[13], dx[14]});
}

// ---------------------------------------------------------------------------
// injectAndResetCovariance(): the complete ESEKF injection step.
//
// Injecting dtheta moves the nominal attitude, so the residual error is now
// parametrized about a DIFFERENT nominal, and the covariance has to be mapped
// across with it: P <- G P G^T.
//
// Derivation for THIS filter's convention (right/local error, q_t = q (x)
// Exp(dtheta)). Write the injected correction as `a` and the true error as
// `a + eps`. After injection the new nominal is q' = q (x) Exp(a), so the
// residual error re-expressed about q' is
//
//   Exp(dtheta_new) = Exp(-a) (x) Exp(a + eps) = Exp( J_r(a) * eps )
//
// giving G_theta = J_r(a) = I - 0.5*[a]x + O(|a|^2), where J_r is the right
// Jacobian of SO(3). The remaining 12 states are additive, so G is
// blkdiag(I - 0.5*[a]x, I_12) and only rows/columns 0..2 of P change.
//
// NOTE the sign. Several references quote G = I + 0.5*[a]x; that is the LEFT
// (global/world) error convention. Using it here would apply the reset
// backwards, doubling the error instead of removing it. Verified numerically
// against a finite-difference derivation of the reset map.
//
// The correction is O(|a|/2) on the attitude block -- a percent-level effect
// for typical per-update corrections, not a stability issue, but it is what
// makes the prediction -> update -> injection -> reset chain a consistent
// ESEKF rather than an EKF with a quaternion glued on.
// ---------------------------------------------------------------------------
void ESEKF::injectAndResetCovariance(const float dx[ESEKF_STATE_DIM])
{
	injectErrorState(dx);

	// G_theta = I - 0.5 * [dtheta]x
	const float hx = 0.5f * dx[0];
	const float hy = 0.5f * dx[1];
	const float hz = 0.5f * dx[2];

	const float G[3][3] = {
		{ 1.0f,  hz,   -hy  },
		{ -hz,   1.0f,  hx  },
		{  hy,  -hx,    1.0f}
	};

	// P <- G P G^T with G = blkdiag(G_theta, I12): only the first three rows
	// and columns are touched, and every output entry depends on exactly three
	// input entries -- the same column in the row pass, the same row in the
	// column pass. Three scalars therefore carry the whole dependency and no
	// scratch array is needed. That matters because this runs nested inside
	// kalmanUpdate()'s frame: the 360 bytes a 3xN plus Nx3 scratch would cost
	// here stack on top of the caller's, at the deepest point in the filter.
	for (int j = 0; j < ESEKF_STATE_DIM; ++j)
	{
		const float p0 = P[0][j], p1 = P[1][j], p2 = P[2][j];
		P[0][j] = G[0][0] * p0 + G[0][1] * p1 + G[0][2] * p2;
		P[1][j] = G[1][0] * p0 + G[1][1] * p1 + G[1][2] * p2;
		P[2][j] = G[2][0] * p0 + G[2][1] * p1 + G[2][2] * p2;
	}

	for (int i = 0; i < ESEKF_STATE_DIM; ++i)
	{
		const float p0 = P[i][0], p1 = P[i][1], p2 = P[i][2];
		P[i][0] = p0 * G[0][0] + p1 * G[0][1] + p2 * G[0][2];
		P[i][1] = p0 * G[1][0] + p1 * G[1][1] + p2 * G[1][2];
		P[i][2] = p0 * G[2][0] + p1 * G[2][1] + p2 * G[2][2];
	}

	// The error state itself is reset to zero implicitly: dx is a local built
	// fresh on every update and discarded here, so nothing carries over.
}

// ---------------------------------------------------------------------------
// symmetrizeCovariance(): P <- 0.5 * (P + P^T).
//
// P is symmetric in exact arithmetic, but F*P*F^T and the Joseph product are
// each computed entry-by-entry in float32, so P[i][j] and P[j][i] accumulate
// their round-off along different summation orders and drift apart. Once P is
// asymmetric the Kalman gain is no longer the gain of any covariance, and the
// asymmetry feeds back through the next FPF^T. Forcing symmetry is 105 adds
// and costs nothing measurable.
//
// This is a numerical hygiene measure, not a correctness patch: it cannot mask
// a wrong Jacobian, because a wrong Jacobian produces a symmetric-but-wrong P.
// ---------------------------------------------------------------------------
void ESEKF::symmetrizeCovariance()
{
	for (int i = 0; i < ESEKF_STATE_DIM; ++i)
	{
		for (int j = i + 1; j < ESEKF_STATE_DIM; ++j)
		{
			const float m = 0.5f * (P[i][j] + P[j][i]);
			P[i][j] = m;
			P[j][i] = m;
		}
	}
}

// ---------------------------------------------------------------------------
// checkFinite(): the filter's numerical health check, run after every
// propagation and every fusion. Returns true while the filter is healthy.
//
// It does three things, in this order:
//
//  1. Screens the nominal state and the covariance diagonal for non-finite or
//     meaningfully negative entries.
//  2. While healthy, keeps a rolling snapshot of the NOMINAL state. That is the
//     rollback point, and it costs 16 float stores per call.
//  3. On a fault, REPAIRS: rolls the nominal state back to that snapshot,
//     re-seeds P to its initial diagonal, counts the event, and carries on.
//
// Step 3 is the part worth explaining, because latching is the obvious design
// and it is wrong in the air. A latched filter refuses every subsequent
// predict() and update(), so the aircraft has no attitude at all -- and it
// cannot get one back, because recovery means re-alignment, re-alignment reads
// gravity off the accelerometer, and that only works on a stationary airframe.
// An aircraft that has just lost its estimate is the opposite of stationary, so
// the filter would sit dead and silent for the remainder of the flight.
//
// Rolling back one step is strictly better: one predict step of attitude error
// is a millisecond of drift, whereas the alternative is no attitude at all. The
// covariance is NOT rolled back -- a covariance that has gone indefinite is not
// worth restoring, and re-seeding it to the initial diagonal is both cheaper
// (no 900-float snapshot on the hot path) and the honest statement of how much
// the filter still knows: the aiding sources re-tighten it within a second.
//
// diverged_ is therefore reached only when there is no rollback point, i.e. the
// fault came out of initialize() itself. That happens on the ground, before
// anything flies, and there refusing is exactly right.
// ---------------------------------------------------------------------------
bool ESEKF::checkFinite()
{
	bool ok = std::isfinite(q.w) && std::isfinite(q.x) && std::isfinite(q.y) && std::isfinite(q.z) &&
	          std::isfinite(velocity_NED_.x) && std::isfinite(velocity_NED_.y) && std::isfinite(velocity_NED_.z) &&
	          std::isfinite(position_NED_.x) && std::isfinite(position_NED_.y) && std::isfinite(position_NED_.z) &&
	          std::isfinite(gyro_bias.x) && std::isfinite(gyro_bias.y) && std::isfinite(gyro_bias.z) &&
	          std::isfinite(accel_bias.x) && std::isfinite(accel_bias.y) && std::isfinite(accel_bias.z);

	// Scanning the diagonal is enough: any NaN reaching an off-diagonal entry
	// reaches the diagonal on the next FPF^T, and this runs every predict().
	//
	// A NEGATIVE diagonal entry needs a distinction that a plain `< 0.0f` test
	// does not make. A variance of -1e-30 on a state whose variance is 1e-9 is
	// float32 round-off in the Joseph product, not divergence -- and latching
	// diverged_ on it would kill an otherwise healthy filter permanently, in
	// flight, for a rounding error. A variance of -5.0 on the same state is
	// real: the covariance is indefinite and nothing downstream of it means
	// anything. So the test is RELATIVE to the size of the covariance: round-off
	// is clamped away (EKF3 and PX4's ekf2 do the same in ConstrainVariances()),
	// and anything meaningfully negative still latches.
	//
	// Clamping to zero rather than to a floor is safe here because every state
	// has a strictly positive Q diagonal, so the next predict() re-inflates it.
	float max_diag = 0.0f;
	for (int i = 0; i < ESEKF_STATE_DIM; ++i)
	{
		if (std::isfinite(P[i][i]) && P[i][i] > max_diag)
		{
			max_diag = P[i][i];
		}
	}
	const float neg_tol = -1.0e-6f * max_diag;

	for (int i = 0; ok && i < ESEKF_STATE_DIM; ++i)
	{
		if (!std::isfinite(P[i][i]))
		{
			ok = false;
		}
		else if (P[i][i] < neg_tol)
		{
			ok = false; // genuinely indefinite
		}
		else if (P[i][i] < 0.0f)
		{
			P[i][i] = 0.0f; // round-off
		}
	}

	if (ok)
	{
		// Healthy: take the rollback snapshot. Nominal state only -- see the
		// banner above for why P is deliberately not part of it.
		last_good_q_           = q;
		last_good_velocity_    = velocity_NED_;
		last_good_position_    = position_NED_;
		last_good_gyro_bias_   = gyro_bias;
		last_good_accel_bias_  = accel_bias;
		have_last_good_        = true;
		return true;
	}

	if (fault_count_ < 0xFFFFFFFFU)
	{
		fault_count_++;
	}

	if (!have_last_good_)
	{
		// Nothing to roll back to: the fault came from alignment itself. Refuse
		// to fly rather than invent a state.
		diverged_ = true;
		return false;
	}

	// Repair. The nominal state goes back one step; the covariance is re-seeded
	// rather than restored, so the filter is honest about having just lost
	// confidence and the aiding sources tighten it again over the next second.
	q             = last_good_q_;
	velocity_NED_ = last_good_velocity_;
	position_NED_ = last_good_position_;
	gyro_bias     = last_good_gyro_bias_;
	accel_bias    = last_good_accel_bias_;
	angular_rate_ = {0.0f, 0.0f, 0.0f}; // stale rate must not reach a rate controller

	resetCovarianceDefaults();
	dt_last_ = 0.0f;

	return false;
}

// ---------------------------------------------------------------------------
// kalmanUpdate(): EKF update for a 3-element measurement.
//   S  = H P H^T + R           (3x3 innovation covariance)
//   NIS = nu^T S^-1 nu         (consistency gate)
//   K  = P H^T S^-1            (15x3 Kalman gain)
//   dx = K * innovation        (15x1, local -- injected then discarded)
//   P  = (I-KH) P (I-KH)^T + K R K^T      <-- Joseph form
//
// Returns false if the update was rejected (singular S, or gated out).
//
// On the covariance form: the previous P = (I-KH)P is algebraically correct
// ONLY at the exact optimal gain, and it is not symmetric by construction, so
// float32 round-off accumulates asymmetrically and can in principle drive P
// indefinite. In this filter that failure was hard to provoke -- Q*dt is added
// at 200 Hz with a strictly positive diagonal on all 15 states, which
// continually re-regularizes P and washes the round-off out (a 120 s run held
// max|P-P^T|/max|P| at 1.4e-6 with no negative eigenvalue, even with R driven
// down to 1e-12). So this is insurance, not a bug fix. It is worth taking
// anyway because the conditions that break the simple form are exactly the ones
// a tuned filter drifts toward: a lower predict:update ratio, a Q diagonal
// someone zeroed on a "known" state, or R far tighter than P. Joseph stays
// PSD for ANY gain, which also makes the gated/skipped paths safe.
//
// Joseph is ~1.7x the flops of (I-KH)P here but LESS stack: two 15x15 scratch
// buffers instead of three.
// ---------------------------------------------------------------------------
bool ESEKF::kalmanUpdate(const float innovation[3], const float H[3][ESEKF_STATE_DIM], const float R[3][3],
                          bool force_fuse)
{
	last_test_ratio_   = 0.0f;
	last_update_gated_ = false;

	if (!initialized_ || diverged_)
	{
		return false;
	}
	if (!std::isfinite(innovation[0]) || !std::isfinite(innovation[1]) || !std::isfinite(innovation[2]))
	{
		return false;
	}

	// PHt = P * H^T  (15x3)
	float PHt[ESEKF_STATE_DIM][3];
	for (int i = 0; i < ESEKF_STATE_DIM; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			float sum = 0.0f;
			for (int k = 0; k < ESEKF_STATE_DIM; ++k)
			{
				sum += P[i][k] * H[j][k]; // H^T[k][j] == H[j][k]
			}
			PHt[i][j] = sum;
		}
	}

	// S = H * PHt + R  (3x3)
	float S[3][3];
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			float sum = 0.0f;
			for (int k = 0; k < ESEKF_STATE_DIM; ++k)
			{
				sum += H[i][k] * PHt[k][j];
			}
			S[i][j] = sum + R[i][j];
		}
	}

	float S_inv[3][3];
	if (!matrixInverse3x3(S, S_inv))
	{
		return false; // innovation covariance singular -- skip rather than corrupt P
	}

	// --- Innovation consistency gate ---
	// NIS = nu^T S^-1 nu. A sample far outside its own predicted uncertainty is
	// far more likely to be a bad measurement (GPS multipath, a magnetic
	// transient, a pressure spike) than a genuine state error, and fusing it
	// puts that error straight into position/attitude with no way back.
	//
	// last_test_ratio_ is NIS normalized so 1.0 sits exactly on the gate, which
	// is what the caller-facing wrappers publish for logging. EKF3 exposes the
	// same normalization as posTestRatio / velTestRatio / hgtTestRatio.
	if (nis_gate_ > 0.0f)
	{
		float nis = 0.0f;
		for (int i = 0; i < 3; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				nis += innovation[i] * S_inv[i][j] * innovation[j];
			}
		}
		last_test_ratio_ = nis / nis_gate_;

		if (!(nis <= nis_gate_)) // written so a NaN NIS also rejects
		{
			last_update_gated_ = true;
			if (!force_fuse || !std::isfinite(nis))
			{
				return false;
			}
			// force_fuse: the caller has established that this source has been
			// gated out for so long that continuing to reject it is the greater
			// risk. Fuse anyway. See ESEKF_GATE_RECOVERY_TIMEOUT.
		}
	}

	// K = PHt * S_inv  (15x3)
	float K[ESEKF_STATE_DIM][3];
	for (int i = 0; i < ESEKF_STATE_DIM; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			float sum = 0.0f;
			for (int k = 0; k < 3; ++k)
			{
				sum += PHt[i][k] * S_inv[k][j];
			}
			K[i][j] = sum;
		}
	}

	// dx = K * innovation  (15x1, local)
	float dx[ESEKF_STATE_DIM];
	for (int i = 0; i < ESEKF_STATE_DIM; ++i)
	{
		dx[i] = K[i][0] * innovation[0] + K[i][1] * innovation[1] + K[i][2] * innovation[2];
	}

	// --- Joseph covariance update: P = A P A^T + K R K^T, with A = I - K H ---
	//
	// Formed WITHOUT ever materializing A or A*P as 15x15 matrices. A is the
	// identity plus a rank-3 correction, so the whole update is a sequence of
	// rank-3 in-place operations on P using only 15x3 scratch:
	//
	//   A P        = P - K (H P)
	//   (A P) A^T  = (A P) - ((A P) H^T) K^T
	//
	// That is exact Joseph -- not an approximation of it -- at 2.9k flops and
	// ~360 bytes of scratch instead of 6.8k flops and 1.8 kB. The dense version
	// pushed this function's frame to ~3.4 kB, which on an STM32 with typical
	// 1-2 kB FreeRTOS task stacks is a stack overflow, not a performance note.

	// One 45-float scratch serving two different shapes in sequence: first
	// HP (3x15), then -- once HP has been fully consumed by the P <- A P pass
	// -- W (15x3). The union makes the reuse explicit and keeps both index
	// expressions natural; each member is completely written before it is
	// read, so no member is ever read while inactive. Two separate buffers
	// would add 180 bytes at the deepest point of the filter's call graph.
	union Scratch45
	{
		float HP[3][ESEKF_STATE_DIM];
		float W[ESEKF_STATE_DIM][3];
	} sc;

	// HP = H * P  (3x15), snapshotted from P before P is touched.
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < ESEKF_STATE_DIM; ++j)
		{
			float sum = 0.0f;
			for (int k = 0; k < ESEKF_STATE_DIM; ++k)
			{
				sum += H[i][k] * P[k][j];
			}
			sc.HP[i][j] = sum;
		}
	}

	// P <- A P = P - K * HP
	for (int i = 0; i < ESEKF_STATE_DIM; ++i)
	{
		for (int j = 0; j < ESEKF_STATE_DIM; ++j)
		{
			P[i][j] -= K[i][0] * sc.HP[0][j] + K[i][1] * sc.HP[1][j] + K[i][2] * sc.HP[2][j];
		}
	}

	// W = (A P) H^T - K R   (15x3). The final rank-3 correction is then a
	// single product, P <- (A P) - W K^T, because
	//   A P A^T + K R K^T = (A P) - ((A P) H^T) K^T + (K R) K^T
	//                     = (A P) - ((A P) H^T - K R) K^T.
	// Folding K*R into W instead of giving it its own 15x3 buffer costs nothing
	// -- the same multiply-accumulates happen -- and removes 180 bytes from the
	// deepest stack frame the filter ever reaches.
	// (HP is dead from here; the same storage becomes W.)
	for (int i = 0; i < ESEKF_STATE_DIM; ++i)
	{
		for (int k = 0; k < 3; ++k)
		{
			float sum = 0.0f;
			for (int j = 0; j < ESEKF_STATE_DIM; ++j)
			{
				sum += P[i][j] * H[k][j]; // H^T[j][k] == H[k][j]
			}
			sum -= K[i][0] * R[0][k] + K[i][1] * R[1][k] + K[i][2] * R[2][k];
			sc.W[i][k] = sum;
		}
	}

	// P <- (A P) - W K^T
	for (int i = 0; i < ESEKF_STATE_DIM; ++i)
	{
		for (int j = 0; j < ESEKF_STATE_DIM; ++j)
		{
			P[i][j] -= sc.W[i][0] * K[j][0] + sc.W[i][1] * K[j][1] + sc.W[i][2] * K[j][2];
		}
	}

	// Inject AFTER the covariance update: P above is the covariance of the
	// error about the pre-injection nominal, and injectAndResetCovariance()
	// then maps it onto the post-injection nominal via P <- G P G^T.
	injectAndResetCovariance(dx);

	symmetrizeCovariance();
	checkFinite();
	return true;
}

// ---------------------------------------------------------------------------
// kalmanUpdateScalar(): Joseph update for a 1-D measurement whose Jacobian has
// a single non-zero entry, z = h * x[state_index].
//
// The barometer previously reused the generic path's structure and built two
// full 15x15 matrices to fuse ONE number: 3375 flops and 1.8 kB of stack for a
// rank-1 update. Exploiting the structure makes it O(N^2) (~500 flops) with a
// 60-byte scratch, which matters when this runs at 20-50 Hz next to everything
// else on the bus.
//
//   S   = h^2 * P[i][i] + r
//   K   = h * P[:,i] / S
//   A   = I - K H  has one non-zero column: A[:,i] = e_i - h*K
//   P   = A P A^T + K r K^T
// ---------------------------------------------------------------------------
bool ESEKF::kalmanUpdateScalar(float innovation, int state_index, float h, float r, bool force_fuse)
{
	last_test_ratio_   = 0.0f;
	last_update_gated_ = false;

	if (!initialized_ || diverged_)
	{
		return false;
	}
	if (state_index < 0 || state_index >= ESEKF_STATE_DIM)
	{
		return false;
	}
	if (!std::isfinite(innovation) || !std::isfinite(h) || !std::isfinite(r))
	{
		return false;
	}

	const int i0 = state_index;

	const float S = h * h * P[i0][i0] + r;
	if (!(S > 1.0e-12f)) // also rejects NaN and a non-positive S
	{
		return false;
	}
	const float S_inv = 1.0f / S;

	if (nis_gate_ > 0.0f)
	{
		const float nis = innovation * innovation * S_inv;
		last_test_ratio_ = nis / nis_gate_;

		if (!(nis <= nis_gate_))
		{
			last_update_gated_ = true;
			if (!force_fuse || !std::isfinite(nis))
			{
				return false;
			}
		}
	}

	float K[ESEKF_STATE_DIM];
	for (int i = 0; i < ESEKF_STATE_DIM; ++i)
	{
		K[i] = h * P[i][i0] * S_inv;
	}

	float dx[ESEKF_STATE_DIM];
	for (int i = 0; i < ESEKF_STATE_DIM; ++i)
	{
		dx[i] = K[i] * innovation;
	}

	// A = I - K H. H is h*e_i0^T, so A = I - h * K * e_i0^T: identity except
	// column i0, which is e_i0 - h*K. Hence
	//   (A P A^T)[i][j] = P[i][j] + a_i * P[i0][j] + a_j * P[i][i0]
	//                     + a_i * a_j * P[i0][i0],   with a_i = -h*K[i].
	float a[ESEKF_STATE_DIM];
	for (int i = 0; i < ESEKF_STATE_DIM; ++i)
	{
		a[i] = -h * K[i];
	}
	const float Pii = P[i0][i0];

	float Pcol[ESEKF_STATE_DIM]; // P[:,i0] before it is overwritten
	for (int i = 0; i < ESEKF_STATE_DIM; ++i)
	{
		Pcol[i] = P[i][i0];
	}

	for (int i = 0; i < ESEKF_STATE_DIM; ++i)
	{
		for (int j = 0; j < ESEKF_STATE_DIM; ++j)
		{
			P[i][j] = P[i][j]
			        + a[i] * Pcol[j]          // P[i0][j] == P[j][i0] == Pcol[j] (P symmetric)
			        + a[j] * Pcol[i]
			        + a[i] * a[j] * Pii
			        + K[i] * r * K[j];        // + K R K^T
		}
	}

	injectAndResetCovariance(dx);

	symmetrizeCovariance();
	checkFinite();
	return true;
}

// ---------------------------------------------------------------------------
// gateRecoveryDue(): has this source been gated out long enough that we must
// fuse the next sample regardless of the innovation gate?
//
// An innovation gate protects the state from a bad measurement, but it is a
// positive feedback loop: reject a sample because the state disagrees with it,
// and the state is then free to drift further, making the next sample look even
// worse. One large transient can therefore lock a perfectly healthy source out
// permanently while the filter dead-reckons beside it.
//
// EKF3 breaks the loop by fusing on `(posCheckPassed || posTimeout ||
// badIMUdata)` -- a timed-out source is fused whether or not it passes. Same
// idea here: once a source has been continuously rejected for longer than
// ESEKF_GATE_RECOVERY_TIMEOUT, trust the sensor over the state.
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// kalmanUpdateScalarVec3(): Joseph update for a scalar measurement whose
// Jacobian is non-zero on one 3-state block only.
//
//   S   = H P H^T + r                       (scalar)
//   K   = P H^T / S                          (15x1)
//   P  <- (I-KH) P (I-KH)^T + K r K^T        (Joseph)
//
// Written out rather than reduced to P - K(HP): with the optimal gain the two
// are algebraically identical, but only the Joseph expansion stays symmetric
// and positive semi-definite under float32 rounding, and this filter has no
// second chance if P goes indefinite in flight.
// ---------------------------------------------------------------------------
bool ESEKF::kalmanUpdateScalarVec3(float innovation, int first_state, const float h3[3],
                                    float r, bool force_fuse)
{
	last_test_ratio_   = 0.0f;
	last_update_gated_ = false;

	if (!initialized_ || diverged_)
	{
		return false;
	}
	if (first_state < 0 || first_state + 3 > ESEKF_STATE_DIM)
	{
		return false;
	}
	if (!std::isfinite(innovation) || !std::isfinite(r) ||
	    !std::isfinite(h3[0]) || !std::isfinite(h3[1]) || !std::isfinite(h3[2]))
	{
		return false;
	}

	const int i0 = first_state;

	// hp = H P  (1x15). H is zero outside the block, so this is three rows.
	float hp[ESEKF_STATE_DIM];
	for (int j = 0; j < ESEKF_STATE_DIM; ++j)
	{
		hp[j] = h3[0] * P[i0][j] + h3[1] * P[i0 + 1][j] + h3[2] * P[i0 + 2][j];
	}

	const float hph = h3[0] * hp[i0] + h3[1] * hp[i0 + 1] + h3[2] * hp[i0 + 2];
	const float S   = hph + r;
	if (!(S > 1.0e-12f)) // also rejects NaN and a non-positive S
	{
		return false;
	}
	const float S_inv = 1.0f / S;

	if (nis_gate_ > 0.0f)
	{
		const float nis = innovation * innovation * S_inv;
		last_test_ratio_ = nis / nis_gate_;

		if (!(nis <= nis_gate_))
		{
			last_update_gated_ = true;
			if (!force_fuse || !std::isfinite(nis))
			{
				return false;
			}
		}
	}

	// K = P H^T / S. (P H^T)_i == (H P)_i because P is symmetric.
	float K[ESEKF_STATE_DIM];
	float dx[ESEKF_STATE_DIM];
	for (int i = 0; i < ESEKF_STATE_DIM; ++i)
	{
		K[i]  = hp[i] * S_inv;
		dx[i] = K[i] * innovation;
	}

	// (A P A^T + K r K^T)[i][j]
	//     = P[i][j] - K_i hp[j] - K_j hp[i] + K_i K_j (hph + r)
	for (int i = 0; i < ESEKF_STATE_DIM; ++i)
	{
		for (int j = 0; j < ESEKF_STATE_DIM; ++j)
		{
			P[i][j] += -K[i] * hp[j] - K[j] * hp[i] + K[i] * K[j] * S;
		}
	}

	injectAndResetCovariance(dx);

	symmetrizeCovariance();
	checkFinite();
	return true;
}

bool ESEKF::sourceTimedOut(double last_fuse_t, float timeout) const
{
	// The difference is taken in double (the clock's type) and only then
	// narrowed: elapsed times are small, so the float result is exact, while
	// the subtraction stays immune to the clock's absolute magnitude.
	return (float)(filter_time_ - last_fuse_t) > timeout;
}

// ---------------------------------------------------------------------------
// resetStateBlockCovariance(): zero rows/cols [first,last] of P, then seed the
// diagonal. This is the covariance half of a state reset -- the old covariance
// described uncertainty about a state we have just thrown away, and its
// cross-correlations to the surviving states are no longer meaningful.
// Mirrors EKF3's zeroStatesVarCov().
// ---------------------------------------------------------------------------
void ESEKF::resetStateBlockCovariance(int first, int last, float variance)
{
	if (first < 0 || last >= ESEKF_STATE_DIM || first > last)
	{
		return;
	}

	// A seeded variance of exactly zero is not "perfect knowledge", it is a
	// dead state: the row and column of P are zero, so K is zero for that
	// state on every subsequent update and no measurement can ever move it
	// again -- and Q only re-inflates the states it has a diagonal entry for.
	// A caller passing 0, a negative, or a NaN (a parameter read back from a
	// corrupt config, say) gets the floor instead of a silently frozen state.
	const float seed = (std::isfinite(variance) && variance > kMinSeedVariance)
	                     ? variance
	                     : kMinSeedVariance;

	for (int i = first; i <= last; ++i)
	{
		for (int j = 0; j < ESEKF_STATE_DIM; ++j)
		{
			P[i][j] = 0.0f;
			P[j][i] = 0.0f;
		}
		P[i][i] = seed;
	}
}

void ESEKF::resetPositionTo(const Vector3f &position_ned, float variance)
{
	if (!std::isfinite(position_ned.x) || !std::isfinite(position_ned.y) || !std::isfinite(position_ned.z))
	{
		return;
	}

	// Record the discontinuity before applying it. A controller integrating
	// position has to shift its own reference by this delta or it will read the
	// jump as a tracking error and fight it.
	last_pos_reset_delta_ = vectorSub(position_ned, position_NED_);
	last_pos_reset_t_     = filter_time_;

	// Classify the jump. EKF3's EK3_GLITCH_RAD is a HORIZONTAL radius, so the
	// down component is deliberately excluded -- a vertical re-centering after
	// a baro timeout is a different event with a different failsafe.
	const float dn = last_pos_reset_delta_.x;
	const float de = last_pos_reset_delta_.y;
	last_pos_reset_exceeded_glitch_ =
		(glitch_radius_ > 0.0f) && ((dn * dn + de * de) > (glitch_radius_ * glitch_radius_));

	position_NED_ = position_ned;
	resetStateBlockCovariance(6, 8, variance);
}

void ESEKF::resetVelocityTo(const Vector3f &velocity_ned, float variance)
{
	if (!std::isfinite(velocity_ned.x) || !std::isfinite(velocity_ned.y) || !std::isfinite(velocity_ned.z))
	{
		return;
	}
	last_vel_reset_delta_ = vectorSub(velocity_ned, velocity_NED_);
	last_vel_reset_t_     = filter_time_;
	velocity_NED_ = velocity_ned;
	resetStateBlockCovariance(3, 5, variance);
}

void ESEKF::resetVerticalPositionTo(float down_m, float variance)
{
	if (!std::isfinite(down_m))
	{
		return;
	}
	last_pos_reset_delta_ = Vector3f{0.0f, 0.0f, down_m - position_NED_.z};
	last_pos_reset_t_     = filter_time_;
	// This path refreshes the reset timestamp, so it must also settle the
	// glitch flag: leaving a `true` from an earlier horizontal glitch would
	// re-arm gps_glitching for another full window on a purely vertical event.
	last_pos_reset_exceeded_glitch_ = false;
	position_NED_.z = down_m;
	resetStateBlockCovariance(8, 8, variance);
}

bool ESEKF::gateRecoveryDue(double last_fuse_t, double last_reject_t) const
{
	if (last_reject_t <= kNeverFused)
	{
		return false; // nothing has been rejected yet
	}
	// Time since the source last actually contributed. A source that has never
	// been fused but is being rejected also qualifies.
	const float since_fuse = (float)(filter_time_ - last_fuse_t);
	return since_fuse > ESEKF_GATE_RECOVERY_TIMEOUT;
}

// ---------------------------------------------------------------------------
// updateAccelerometer(): treats gravity as a weak "up vector" reference to
// correct roll/pitch drift. Gated on |accel| being close to local gravity
// magnitude -- when the vehicle is under real acceleration (a maneuver, not
// just gravity), the reading is not a trustworthy "up" reference and the
// update is skipped rather than corrupting attitude/accel_bias with it.
// Tune the gate width with setAccelGateThreshold().
// Predicted measurement: z_pred = R^T * (-g) + accel_bias
// Jacobian: d(z)/d(dtheta) = skew(v), v = R^T * (-g); d(z)/d(dba) = I
// ---------------------------------------------------------------------------
bool ESEKF::updateAccelerometer(const Vector3f &accel)
{
	if (!initialized_ || diverged_)
	{
		return false;
	}
	if (!std::isfinite(accel.x) || !std::isfinite(accel.y) || !std::isfinite(accel.z))
	{
		return false;
	}

	// How far this reading is from being pure gravity. Recorded before the gate
	// so getAccelMotion() reports the vibration floor even while every sample is
	// being rejected -- which is exactly the case the number has to diagnose.
	const float accel_norm   = vectorNorm(accel);
	const float gravity_norm = vectorNorm(g);
	accel_motion_ = std::fabs(accel_norm - gravity_norm);

	// The gate draws the line where the reading stops being a gravity reference
	// at all. Inside it, trust is graded rather than binary: R is inflated with
	// the square of how far into the gate the sample sits, so a marginal sample
	// is de-weighted instead of being taken at full confidence one cycle and
	// thrown away the next. A bare threshold makes the filter chatter between
	// full aiding and none for a vehicle hovering right at the limit, and that
	// shows up as attitude that twitches with throttle.
	float R_use[3][3];
	if (accel_gate_threshold_ > 0.0f)
	{
		if (accel_motion_ > accel_gate_threshold_)
		{
			return false; // real acceleration -- skip, don't corrupt attitude/bias
		}

		const float ratio = accel_motion_ / accel_gate_threshold_; // [0, 1]
		accel_inflation_ = 1.0f + kAccelInflationGain * ratio * ratio;
	}
	else
	{
		accel_inflation_ = 1.0f; // gate disabled: no grading either
	}

	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			R_use[i][j] = R_accel[i][j] * accel_inflation_;
		}
	}

	float R[3][3];
	quaternionToRotationMatrix(q, R);

	// v = R^T * (-g)
	Vector3f neg_g{-g.x, -g.y, -g.z};
	Vector3f v{
		R[0][0] * neg_g.x + R[1][0] * neg_g.y + R[2][0] * neg_g.z,
		R[0][1] * neg_g.x + R[1][1] * neg_g.y + R[2][1] * neg_g.z,
		R[0][2] * neg_g.x + R[1][2] * neg_g.y + R[2][2] * neg_g.z
	};

	Vector3f z_pred = vectorAdd(v, accel_bias);
	Vector3f innov  = vectorSub(accel, z_pred);
	float innovation[3] = {innov.x, innov.y, innov.z};

	float S_v[3][3];
	skewSymmetric(v, S_v);

	float H[3][ESEKF_STATE_DIM];
	std::memset(H, 0, sizeof(H));
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			H[i][j] = S_v[i][j];
		}
		H[i][12 + i] = 1.0f; // d(z)/d(dba)
	}

	return kalmanUpdate(innovation, H, R_use);
}

// ---------------------------------------------------------------------------
// Filter status / integrity reporting.
//
// This is the interface an autonomous flight controller actually flies on. The
// numbers above it -- position, velocity, attitude -- are only meaningful when
// read together with these flags, because the estimator will keep producing
// smooth, plausible, confidently-wrong output long after it has stopped being
// aided by anything. That is not a defect to be fixed; it is what an inertial
// filter DOES when its absolute references go away, and the only defence is to
// say so out loud.
//
// Shaped after EKF3's nav_filter_status, and like EKF3 every flag is driven by
// the time since a source was last SUCCESSFULLY FUSED -- never by whether data
// is arriving. A receiver streaming fixes that the innovation gate rejects on
// every sample is contributing exactly nothing, and a health check that watches
// "GPS connected" or "satellite count" will report everything is fine.
// ---------------------------------------------------------------------------
ESEKFStatus ESEKF::getStatus() const
{
	ESEKFStatus st;

	const float since_gps_pos = (float)(filter_time_ - last_gps_pos_fuse_t_);
	const float since_gps_vel = (float)(filter_time_ - last_gps_vel_fuse_t_);
	const float since_baro    = (float)(filter_time_ - last_baro_fuse_t_);
	const float since_mag     = (float)(filter_time_ - last_mag_fuse_t_);

	const bool healthy = initialized_ && !diverged_;

	// A GPS position or velocity solution fused inside the aiding window is
	// what makes horizontal position absolute.
	const bool gps_pos_aiding = (since_gps_pos < ESEKF_GPS_AID_TIMEOUT);
	const bool gps_vel_aiding = (since_gps_vel < ESEKF_GPS_AID_TIMEOUT);
	const bool hgt_aiding     = (since_baro < ESEKF_HGT_AID_TIMEOUT) || gps_pos_aiding;

	st.attitude_valid  = healthy;
	st.horiz_vel_valid = healthy && (gps_vel_aiding || gps_pos_aiding);
	st.vert_vel_valid  = healthy && (gps_vel_aiding || hgt_aiding);
	st.horiz_pos_valid = healthy && gps_pos_aiding;
	st.vert_pos_valid  = healthy && hgt_aiding;

	// using_gps uses the shorter window: EKF3 reports it on
	// (now - lastGpsPosPassTime_ms) < 4000.
	st.using_gps = healthy && (since_gps_pos < ESEKF_GPS_USE_TIMEOUT);

	// DEAD RECKONING: no absolute horizontal position source is being fused.
	// Position is now pure inertial integration and its error grows without
	// bound -- roughly quadratically, driven by residual tilt and accel bias.
	// A vehicle in this state should be heading for a failsafe, not a waypoint.
	st.dead_reckoning = healthy && !gps_pos_aiding && !gps_vel_aiding;

	// GPS glitching: the most recent GPS activity was a REJECTION, and it was
	// recent. Deliberately not conditioned on using_gps -- the first seconds of
	// a glitch are exactly when using_gps is still true (the last good fix is
	// only moments old), which is precisely the window a flight controller needs
	// warning in.
	//
	// A position reset larger than the glitch radius counts as a glitch too.
	// That is the whole point of setGlitchRadius(): a routine few-metre
	// re-centering after a brief outage and a 500 m teleport both raise
	// position_reset, but only the second one means the vehicle's idea of
	// where it is just changed by more than it can safely absorb, and only the
	// second one should trip a failsafe.
	const float since_gps_reject = (float)(filter_time_ - last_gps_pos_reject_t_);
	const float since_pos_reset  = (float)(filter_time_ - last_pos_reset_t_);
	const bool  glitch_reset     = last_pos_reset_exceeded_glitch_ &&
	                               (last_pos_reset_t_ > kNeverFused) &&
	                               (since_pos_reset < ESEKF_GPS_USE_TIMEOUT);

	st.gps_glitching = healthy &&
	                    (((since_gps_reject < ESEKF_GPS_USE_TIMEOUT) &&
	                      (last_gps_pos_reject_t_ > last_gps_pos_fuse_t_)) ||
	                     glitch_reset);

	// A reset means the position estimate TELEPORTED. Reported for the full
	// aiding window so a slow control loop cannot miss it.
	st.position_reset = healthy &&
	                     (last_pos_reset_t_ > kNeverFused) &&
	                     (since_pos_reset < ESEKF_GPS_USE_TIMEOUT);

	st.mag_aiding  = healthy && (since_mag < ESEKF_MAG_AID_TIMEOUT);
	st.initialized = initialized_;
	st.diverged    = diverged_;

	return st;
}

bool ESEKF::isDeadReckoning() const { return getStatus().dead_reckoning; }

float ESEKF::getTimeSinceGPSFusion() const  { return (float)(filter_time_ - last_gps_pos_fuse_t_); }
float ESEKF::getTimeSinceBaroFusion() const { return (float)(filter_time_ - last_baro_fuse_t_); }
float ESEKF::getTimeSinceMagFusion() const  { return (float)(filter_time_ - last_mag_fuse_t_); }

float ESEKF::getGPSPosTestRatio() const { return pos_test_ratio_; }
float ESEKF::getGPSVelTestRatio() const { return vel_test_ratio_; }
float ESEKF::getBaroTestRatio() const   { return hgt_test_ratio_; }
float ESEKF::getMagTestRatio() const    { return mag_test_ratio_; }

// 1-sigma horizontal position uncertainty, sqrt(P_NN + P_EE). A navigation mode
// should refuse to engage, and an engaged one should degrade, above a threshold.
float ESEKF::getHorizontalPositionError() const
{
	const float v = P[6][6] + P[7][7];
	return (v > 0.0f) ? std::sqrt(v) : 0.0f;
}

float ESEKF::getVerticalPositionError() const
{
	return (P[8][8] > 0.0f) ? std::sqrt(P[8][8]) : 0.0f;
}

Vector3f ESEKF::getLastPositionResetDelta() const { return last_pos_reset_delta_; }
Vector3f ESEKF::getLastVelocityResetDelta() const { return last_vel_reset_delta_; }
float ESEKF::getLastPositionResetTime() const     { return (float)last_pos_reset_t_; }
float ESEKF::getLastVelocityResetTime() const     { return (float)last_vel_reset_t_; }

// A non-finite or negative radius would make every comparison against it
// false, quietly disabling glitch reporting; clamp to 0 (= report nothing)
// only when that is what the caller actually asked for.
void ESEKF::setGlitchRadius(float radius_m)
{
	glitch_radius_ = (std::isfinite(radius_m) && radius_m > 0.0f) ? radius_m : 0.0f;
}
float ESEKF::getGlitchRadius() const        { return glitch_radius_; }

// ---------------------------------------------------------------------------
// updateMagnetometer(): corrects yaw (and, weakly, roll/pitch) using the
// reference field mag_ref_ computed at inilialize(). No mag bias state is
// modeled here.
// Predicted measurement: z_pred = R^T * mag_ref_
// Jacobian: d(z)/d(dtheta) = skew(v), v = R^T * mag_ref_
// ---------------------------------------------------------------------------
bool ESEKF::updateMagnetometer(const Vector3f &mag)
{
	if (!initialized_ || diverged_)
	{
		return false;
	}
	if (!std::isfinite(mag.x) || !std::isfinite(mag.y) || !std::isfinite(mag.z))
	{
		return false;
	}

	// ---- HEADING ONLY. The magnetometer must never move roll or pitch. ----
	//
	// This used to fuse the field as a 3-vector, z = R^T * mag_ref_, with
	// H = [mag_ref_ in body]x on the attitude block. That is the textbook
	// construction and it is WRONG for a filter shaped like this one, because
	// mag_ref_ is fixed at alignment and there are no earth-field states to
	// absorb the difference between it and the field actually present. Every
	// discrepancy -- a local anomaly, a steel gate, current through a nearby
	// harness, flying somewhere with a different inclination than where the
	// vehicle was switched on -- was therefore fused as a full 3-axis attitude
	// error, and two thirds of that correction lands in ROLL AND PITCH.
	//
	// A vehicle told it is banked when it is level does exactly the wrong thing.
	// The accelerometer fights the error back, so the symptom is not a clean
	// failure: it is attitude that leans toward magnetic clutter and recovers
	// slowly, which is very hard to attribute in flight.
	//
	// Geomagnetism carries no information about tilt anyway -- gravity already
	// fixes two axes, and the only thing the magnetometer is needed for is the
	// third. So fuse exactly that: one scalar, the heading discrepancy, with a
	// Jacobian that can only rotate the estimate about the EARTH VERTICAL. Roll
	// and pitch are then structurally untouchable by the magnetometer, whatever
	// it reads. This is what ArduPilot's fuseEulerYaw() does; the 3-axis
	// alternative is only sound alongside the 6 earth/body field states EKF3
	// carries to make the reference observable, which is out of scope here.
	const float mag_norm = vectorNorm(mag);
	if (mag_norm < 1.0e-9f)
	{
		return false; // no field, no heading
	}

	float R[3][3];
	quaternionToRotationMatrix(q, R);

	// Measured field rotated into NED through the CURRENT attitude estimate. If
	// the estimate's yaw is wrong by eps, this lands rotated by eps about down
	// from mag_ref_ -- which is precisely the quantity to null.
	const Vector3f m_ned{
		R[0][0] * mag.x + R[0][1] * mag.y + R[0][2] * mag.z,
		R[1][0] * mag.x + R[1][1] * mag.y + R[1][2] * mag.z,
		R[2][0] * mag.x + R[2][1] * mag.y + R[2][2] * mag.z};

	// Horizontal components only: the vertical one is inclination, which is what
	// must not be allowed to influence the answer.
	const float h_meas = std::sqrt(m_ned.x * m_ned.x + m_ned.y * m_ned.y);
	const float h_ref  = std::sqrt(mag_ref_.x * mag_ref_.x + mag_ref_.y * mag_ref_.y);

	// Near-vertical field (magnetic poles, or a reference set with no horizontal
	// component) carries no heading. Refusing keeps a meaningless angle out of
	// the state; the aiding timeout then reports yaw as unaided, which is true.
	const float kMinHorizField = 1.0e-4f; // fraction-of-Earth-field scale guard
	if (h_meas < kMinHorizField * mag_norm || h_ref < 1.0e-9f)
	{
		return false;
	}

	// eps > 0 means the estimate's yaw is LARGER than the truth.
	float eps = std::atan2(m_ned.y, m_ned.x) - std::atan2(mag_ref_.y, mag_ref_.x);
	eps = wrapPiF(eps);

	// Observation: drive eps to zero. z = 0, h(x) = eps, so innovation = -eps.
	const float innovation = -eps;

	// The attitude error state is a BODY-frame rotation vector (the nominal is
	// updated by right multiplication, see injectErrorState()). A rotation about
	// the earth vertical is therefore u = R^T * [0,0,1], the third ROW of R, and
	// d(eps)/d(dtheta) = u^T. Because u is a unit vector this Jacobian has unit
	// gain in the yaw direction and is exactly zero in every direction
	// orthogonal to it -- the structural guarantee that roll and pitch cannot
	// move, which no amount of tuning could give the 3-axis form.
	const float h3[3] = {R[2][0], R[2][1], R[2][2]};

	// Heading noise, in radians^2, derived from the caller's field-unit noise so
	// setMagNoise() keeps meaning what it says. A field error of sigma across
	// the horizontal component subtends sigma/|H_horizontal| radians of heading,
	// so the angular variance is the field variance divided by H^2. This is why
	// a weak horizontal field -- high latitude, or a badly sited sensor -- is
	// automatically trusted less rather than needing a separate parameter.
	const float mag_var = (R_mag[0][0] + R_mag[1][1] + R_mag[2][2]) * (1.0f / 3.0f);
	float r_yaw = mag_var / (h_meas * h_meas);
	if (!(r_yaw > kMinYawVariance))
	{
		r_yaw = kMinYawVariance; // never claim a perfect heading
	}

	// The magnetometer keeps the forced-fusion escape rather than a state reset.
	// Its innovation is an angle and so is bounded by pi, meaning a mag lockout
	// cannot produce the unbounded error that a position lockout can, and
	// repeatedly fusing is enough to walk yaw back. (EKF3 does have a yaw reset,
	// driven by its GSF yaw estimator, which is well beyond this filter's scope.)
	const bool force = gateRecoveryDue(last_mag_fuse_t_, last_mag_reject_t_);
	const bool fused = kalmanUpdateScalarVec3(innovation, 0, h3, r_yaw, force);

	mag_test_ratio_ = last_test_ratio_;
	if (fused)
	{
		last_mag_fuse_t_ = filter_time_;
	}
	if (last_update_gated_)
	{
		last_mag_reject_t_ = filter_time_;
	}
	return fused;
}

// ---------------------------------------------------------------------------
// updateBarometer(): scalar update on the down-position, z = -p_z.
//
// REFERENCE DATUM -- the important part.
// The residual is formed against baro_ref_, the barometer's OWN reading at
// initialize(), not against gps_ref_.z. That makes this a RELATIVE altitude
// measurement: (altitude - baro_ref_) is the climb above the initialization
// point, and -position_NED_.z is exactly the same quantity, so whatever datum
// the barometer uses cancels out of the residual entirely.
//
// The previous code used gps_ref_.z as the reference, which silently assumes
// the barometer and the GPS share an altitude datum. They do not: barometric
// altitude (pressure altitude against 1013.25 hPa, or QNH) and GPS altitude
// (WGS84 ellipsoidal height, or geoid height) routinely differ by tens to
// hundreds of metres, and that entire offset lands in the residual on the very
// first update. With a 40 m offset this produced a 20 m jump on the first
// barometer sample and settled to a 39.8 m vertical error within 10 s at 20 Hz,
// while GPS pulled the other way and the standoff bled into the vertical
// velocity and accel-bias states through their cross-covariances.
//
// Note the barometer is still treated as UNBIASED about that datum. Real
// barometric altitude drifts with weather over tens of minutes; for long
// flights add a baro-bias state and let GPS observe it, rather than widening
// R_baro.
// ---------------------------------------------------------------------------
bool ESEKF::updateBarometer(float altitude)
{
	if (!initialized_ || diverged_)
	{
		return false;
	}
	if (!std::isfinite(altitude))
	{
		return false;
	}

	const int pz_idx = 8; // index of dp_z within the 15-element error state

	const float rel_alt    = altitude - baro_ref_;  // climb above the init point (up-positive)
	const float z_pred     = -position_NED_.z;      // predicted climb (up-positive)
	const float innovation = rel_alt - z_pred;

	altitude_ = altitude;

	// H is 1x15 with H[pz_idx] = -1 (z = -p_z), zero elsewhere.
	const bool fused = kalmanUpdateScalar(innovation, pz_idx, -1.0f, R_baro, false);

	hgt_test_ratio_ = last_test_ratio_;
	if (fused)
	{
		last_baro_fuse_t_ = filter_time_;
		return true;
	}

	if (last_update_gated_)
	{
		last_baro_reject_t_ = filter_time_;

		if (sourceTimedOut(last_baro_fuse_t_, ESEKF_HGT_AID_TIMEOUT))
		{
			resetVerticalPositionTo(-rel_alt, (R_baro > 0.0f) ? R_baro : 1.0f);
			last_baro_fuse_t_ = filter_time_;
			hgt_test_ratio_   = 0.0f;
			return true;
		}
	}
	return false;
}

// ---------------------------------------------------------------------------
// setGPSReferencePrecise(): sets the local NED origin.
// ---------------------------------------------------------------------------
void ESEKF::setGPSReferencePrecise(double latitude_rad, double longitude_rad, float altitude_m)
{
	// Rejected rather than clamped: an out-of-range latitude means the caller's
	// units are wrong, and clamping would plant the local frame's origin at a
	// pole while every getter kept returning confident numbers. Keeping the
	// previous origin makes the mistake show up as gated GPS and a
	// dead_reckoning status, which is at least visible.
	if (!isValidLatLonRad(latitude_rad, longitude_rad) || !std::isfinite(altitude_m))
	{
		return;
	}

	gps_ref_lat_rad_ = latitude_rad;
	gps_ref_lon_rad_ = longitude_rad;
	gps_ref_ = Vector3f{(float)latitude_rad, (float)longitude_rad, altitude_m};
}

// ---------------------------------------------------------------------------
// updateGPSPrecise(): converts a raw GPS fix to a local NED offset from the
// origin and fuses it as a direct position measurement.
//
// UNITS: latitude/longitude in RADIANS, altitude in metres. Nothing in the
// filter can detect degrees being passed instead -- the result is simply a
// ~57x scale error on the horizontal position -- so the degrees entry point is
// a separate function rather than a runtime guess.
//
// Two things the previous flat-Earth conversion got wrong:
//
//  1. It used the equatorial radius for BOTH axes. North displacement needs the
//     meridional radius of curvature M, which is 0.67% smaller than `a` at the
//     equator and 0.17% smaller at 45 deg -- a pure scale error on the local
//     frame's north axis, ~1.7 m per km travelled at mid latitudes.
//
//  2. It differenced the angles in float32. One float32 ULP at 45 deg latitude
//     (0.785 rad) is 6.0e-8 rad = 0.38 m on the ground, and 0.76 m at 60 deg,
//     so the position measurement arrived pre-quantized more coarsely than a
//     decent receiver's own noise. The differencing is done in double here and
//     only the small local result is narrowed to float, which is exact.
//
// Still a local-tangent-plane approximation: fine to a few km, and the
// linearization point is the origin, so accuracy degrades with distance from it.
// ---------------------------------------------------------------------------
bool ESEKF::updateGPSPrecise(double latitude_rad, double longitude_rad, float altitude_m)
{
	if (!initialized_ || diverged_)
	{
		return false;
	}
	if (!std::isfinite(latitude_rad) || !std::isfinite(longitude_rad) || !std::isfinite(altitude_m))
	{
		return false;
	}

	double r_north, r_east;
	earthRadii(gps_ref_lat_rad_, r_north, r_east);

	const double dlat = latitude_rad - gps_ref_lat_rad_;
	const double dlon = wrapPi(longitude_rad - gps_ref_lon_rad_); // antimeridian-safe
	const float  dalt = altitude_m - gps_ref_.z;

	Vector3f gps_pos_ned{
		(float)(dlat * r_north),   // north
		(float)(dlon * r_east),    // east  (r_east already carries cos(lat))
		-dalt                      // down  (GPS altitude is up-positive, NED z is down)
	};

	Vector3f innov = vectorSub(gps_pos_ned, position_NED_);
	float innovation[3] = {innov.x, innov.y, innov.z};

	float H[3][ESEKF_STATE_DIM];
	std::memset(H, 0, sizeof(H));
	for (int i = 0; i < 3; ++i)
	{
		H[i][6 + i] = 1.0f; // d(position)/d(dp) = I
	}

	const bool fused = kalmanUpdate(innovation, H, R_gps, false);

	pos_test_ratio_ = last_test_ratio_;
	if (fused)
	{
		last_gps_pos_fuse_t_ = filter_time_;
		return true;
	}

	if (last_update_gated_)
	{
		last_gps_pos_reject_t_ = filter_time_;

		// The gate has been rejecting this source for longer than the aiding
		// timeout: the state, not the receiver, is the thing that is wrong.
		// Snap position to the measurement and re-seed its covariance from
		// R_gps rather than trying to fuse the difference away. See
		// resetPositionTo() for why forced fusion is not enough here.
		if (sourceTimedOut(last_gps_pos_fuse_t_, ESEKF_GPS_AID_TIMEOUT))
		{
			const float var = (R_gps[0][0] > 0.0f) ? R_gps[0][0] : 1.0f;
			resetPositionTo(gps_pos_ned, var);
			last_gps_pos_fuse_t_ = filter_time_;
			pos_test_ratio_      = 0.0f;
			return true;
		}
	}
	return false;
}

// ---------------------------------------------------------------------------
// updateGPSVelocity(): fuse the receiver's NED Doppler velocity solution.
//
//   z = v_ned,   H = [0  I  0  0  0],   dz/d(dv) = I
//
// The measurement model is the identity on the velocity states, so the Jacobian
// is exact and needs no linearization -- there is nothing to get wrong here.
//
// Why this matters more than its simplicity suggests: GPS Doppler velocity is
// an independent, direct observation of the velocity states, accurate to
// ~0.05-0.1 m/s where differentiated position is accurate to metres. Because
// the velocity error dynamics couple to attitude through -R[f]x and to accel
// bias through -R (see predictCovariance), an accurate velocity
// measurement is what actually makes TILT and ACCELEROMETER BIAS observable.
// Without it those were the two worst-estimated states in this filter --
// measured accel bias converged to [+0.068, -0.145, +0.013] against a truth of
// [+0.060, -0.090, +0.050]. EKF3 treats GPS velocity as a primary source for
// exactly this reason.
// ---------------------------------------------------------------------------
bool ESEKF::updateGPSVelocity(const Vector3f &velocity_ned)
{
	if (!initialized_ || diverged_)
	{
		return false;
	}
	if (!std::isfinite(velocity_ned.x) || !std::isfinite(velocity_ned.y) || !std::isfinite(velocity_ned.z))
	{
		return false;
	}

	Vector3f innov = vectorSub(velocity_ned, velocity_NED_);
	float innovation[3] = {innov.x, innov.y, innov.z};

	float H[3][ESEKF_STATE_DIM];
	std::memset(H, 0, sizeof(H));
	for (int i = 0; i < 3; ++i)
	{
		H[i][3 + i] = 1.0f; // d(velocity)/d(dv) = I
	}

	const bool fused = kalmanUpdate(innovation, H, R_gps_vel, false);

	vel_test_ratio_ = last_test_ratio_;
	if (fused)
	{
		last_gps_vel_fuse_t_ = filter_time_;
		return true;
	}

	if (last_update_gated_)
	{
		last_gps_vel_reject_t_ = filter_time_;

		if (sourceTimedOut(last_gps_vel_fuse_t_, ESEKF_GPS_AID_TIMEOUT))
		{
			const float var = (R_gps_vel[0][0] > 0.0f) ? R_gps_vel[0][0] : 1.0f;
			resetVelocityTo(velocity_ned, var);
			last_gps_vel_fuse_t_ = filter_time_;
			vel_test_ratio_      = 0.0f;
			return true;
		}
	}
	return false;
}

// gps = {latitude(rad), longitude(rad), altitude(m)}. See updateGPSPrecise()
// for why passing lat/lon through a float32 Vector3f costs ~0.4 m of resolution.
bool ESEKF::updateGPS(const Vector3f &gps)
{
	return updateGPSPrecise((double)gps.x, (double)gps.y, gps.z);
}

// gps_deg = {latitude(deg), longitude(deg), altitude(m)} -- the format almost
// every NMEA / u-blox driver actually produces.
bool ESEKF::updateGPSDegrees(const Vector3f &gps_deg)
{
	const double kDegToRad = 0.017453292519943295;
	return updateGPSPrecise((double)gps_deg.x * kDegToRad,
	                        (double)gps_deg.y * kDegToRad,
	                        gps_deg.z);
}

// ---------------------------------------------------------------------------
// Lifecycle
// ---------------------------------------------------------------------------
void ESEKF::reset()
{
	q = {1.0f, 0.0f, 0.0f, 0.0f};
	velocity_NED_ = {0.0f, 0.0f, 0.0f};
	position_NED_ = {0.0f, 0.0f, 0.0f};
	gyro_bias     = {0.0f, 0.0f, 0.0f};
	accel_bias    = {0.0f, 0.0f, 0.0f};

	// angular_rate_ is filter state, not tuning: leaving the last pre-reset
	// value in place would hand a rate controller a stale body rate on the
	// first cycle after a reset, before any predict() has refreshed it.
	angular_rate_ = {0.0f, 0.0f, 0.0f};

	resetCovarianceDefaults();

	dt_last_     = 0.0f;
	altitude_    = 0.0f;
	baro_ref_    = 0.0f;
	initialized_ = false;
	diverged_    = false;

	// The rollback point describes a state this filter no longer holds.
	// fault_count_ deliberately survives -- it is a record of what happened to
	// this airframe, not filter state, and clearing it would hide a fault that
	// was followed by a reset.
	have_last_good_ = false;

	accel_motion_    = 0.0f;
	accel_inflation_ = 1.0f;

	filter_time_           = 0.0;
	last_gps_pos_fuse_t_   = kNeverFused;
	last_gps_vel_fuse_t_   = kNeverFused;
	last_baro_fuse_t_      = kNeverFused;
	last_mag_fuse_t_       = kNeverFused;
	last_gps_pos_reject_t_ = kNeverFused;
	last_gps_vel_reject_t_ = kNeverFused;
	last_baro_reject_t_    = kNeverFused;
	last_mag_reject_t_     = kNeverFused;
	pos_test_ratio_ = 0.0f;
	vel_test_ratio_ = 0.0f;
	hgt_test_ratio_ = 0.0f;
	mag_test_ratio_ = 0.0f;

	last_pos_reset_delta_ = {0.0f, 0.0f, 0.0f};
	last_vel_reset_delta_ = {0.0f, 0.0f, 0.0f};
	last_pos_reset_t_     = kNeverFused;
	last_vel_reset_t_     = kNeverFused;
	last_pos_reset_exceeded_glitch_ = false;
}

bool ESEKF::isInitialized() const
{
	return initialized_;
}

bool ESEKF::hasDiverged() const
{
	return diverged_;
}

uint32_t ESEKF::getFaultCount() const
{
	return fault_count_;
}

// ---------------------------------------------------------------------------
// State getters / setters
// ---------------------------------------------------------------------------
Quaternionf ESEKF::getOrientation() const { return q; }
Vector3f ESEKF::getEulerAngles() const { return eulerFromQuaternion(q); }
Vector3f ESEKF::getVelocity() const { return velocity_NED_; }
Vector3f ESEKF::getPosition() const { return position_NED_; }

// Inverse of updateGPS()'s flat-Earth conversion: turns the filter's local
// NED position back into {latitude, longitude, altitude} using getGPSReference()
// as the origin. Same unit assumption as updateGPS() (lat/lon in radians) --
// keep the two in sync if that assumption changes.
Vector3f ESEKF::getGPSPosition() const
{
	double r_north, r_east;
	earthRadii(gps_ref_lat_rad_, r_north, r_east);

	const double dlat = (double)position_NED_.x / r_north;  // north -> dlat
	// r_east -> 0 at the poles; clamp so this returns a finite (if meaningless)
	// longitude instead of an Inf that then poisons whatever consumes it.
	const double dlon = (std::fabs(r_east) > 1.0)
	                      ? (double)position_NED_.y / r_east
	                      : 0.0;
	const float  dalt = -position_NED_.z;                   // down -> dalt

	return Vector3f{(float)(gps_ref_lat_rad_ + dlat),
	                (float)(gps_ref_lon_rad_ + dlon),
	                gps_ref_.z + dalt};
}

Vector3f ESEKF::getGyroBias() const { return gyro_bias; }
Vector3f ESEKF::getAccelBias() const { return accel_bias; }
Vector3f ESEKF::getAngularRate() const { return angular_rate_; }
// Derived from position_NED_, not the stale altitude_ member -- predict()
// continuously propagates position but does not touch altitude_, so this
// keeps getAltitude() consistent with getPosition() between barometer updates.
float ESEKF::getAltitude() const { return gps_ref_.z - position_NED_.z; } // absolute altitude = reference + relative climb

// These write the nominal state directly, bypassing every estimator path, so
// they are the shortest route from a bad caller value to a diverged filter.
// A non-finite argument is ignored: quaternionNormalize() cannot rescue a NaN
// quaternion either (the norm is NaN, so its small-norm guard never fires).
void ESEKF::setOrientation(const Quaternionf &q0)
{
	if (!std::isfinite(q0.w) || !std::isfinite(q0.x) ||
	    !std::isfinite(q0.y) || !std::isfinite(q0.z))
	{
		return;
	}
	q = quaternionNormalize(q0);
}

void ESEKF::setVelocity(const Vector3f &v0) { if (isFiniteVec(v0)) velocity_NED_ = v0; }
void ESEKF::setPosition(const Vector3f &p0) { if (isFiniteVec(p0)) position_NED_ = p0; }
void ESEKF::setGyroBias(const Vector3f &b0) { if (isFiniteVec(b0)) gyro_bias = b0; }
void ESEKF::setAccelBias(const Vector3f &b0) { if (isFiniteVec(b0)) accel_bias = b0; }

// ---------------------------------------------------------------------------
// Covariance getters / setters
// ---------------------------------------------------------------------------
void ESEKF::getCovariance(float P_out[ESEKF_STATE_DIM][ESEKF_STATE_DIM]) const
{
	std::memcpy(P_out, P, sizeof(P));
}

// Accepted only if every entry is finite and the diagonal is non-negative --
// the cheap necessary conditions for a covariance. (Full positive-definiteness
// would need an eigen or Cholesky check; that is the caller's business, and
// the Joseph-form updates below stay PSD from any PSD starting point.)
//
// The copy is symmetrized on the way in for the same reason predict() and the
// updates symmetrize on the way out: an asymmetric P is not the covariance of
// anything, and the asymmetry compounds through every subsequent F P F^T.
void ESEKF::setCovariance(const float P_in[ESEKF_STATE_DIM][ESEKF_STATE_DIM])
{
	for (int i = 0; i < ESEKF_STATE_DIM; ++i)
	{
		if (!std::isfinite(P_in[i][i]) || P_in[i][i] < 0.0f)
		{
			return;
		}
		for (int j = 0; j < ESEKF_STATE_DIM; ++j)
		{
			if (!std::isfinite(P_in[i][j]))
			{
				return;
			}
		}
	}

	std::memcpy(P, P_in, sizeof(P));
	symmetrizeCovariance();
}

float ESEKF::getStateVariance(int index) const
{
	if (index < 0 || index >= ESEKF_STATE_DIM)
	{
		return 0.0f;
	}
	return P[index][index];
}

// ---------------------------------------------------------------------------
// Process noise getters / setters
// ---------------------------------------------------------------------------
// Each of these rejects a non-finite or negative sigma^2 outright and keeps
// the previous value. A negative process noise makes P shrink every predict()
// until it goes indefinite; a NaN turns all of P into NaN on the first
// propagation. Neither is recoverable in flight.
void ESEKF::setProcessNoiseGyro(float sigma2)
{
	if (!isValidVariance(sigma2)) return;
	for (int i = 0; i < 3; ++i) Q[i][i] = sigma2;
}

void ESEKF::setProcessNoiseAccel(float sigma2)
{
	if (!isValidVariance(sigma2)) return;
	for (int i = 0; i < 3; ++i) Q[3 + i][3 + i] = sigma2;
}

void ESEKF::setProcessNoiseGyroBias(float sigma2)
{
	if (!isValidVariance(sigma2)) return;
	for (int i = 0; i < 3; ++i) Q[9 + i][9 + i] = sigma2;
}

void ESEKF::setProcessNoiseAccelBias(float sigma2)
{
	if (!isValidVariance(sigma2)) return;
	for (int i = 0; i < 3; ++i) Q[12 + i][12 + i] = sigma2;
}

void ESEKF::getProcessNoise(float Q_out[ESEKF_STATE_DIM][ESEKF_STATE_DIM]) const
{
	std::memcpy(Q_out, Q, sizeof(Q));
}

// ---------------------------------------------------------------------------
// setImuNoiseParameters(): derives Q's attitude/velocity/gyro-bias/accel-bias
// diagonal blocks from IMU datasheet noise-density figures. Units:
//   gyro_noise_density_    : (rad/s)/sqrt(Hz)
//   accel_noise_density_   : (m/s^2)/sqrt(Hz)
//   gyro_bias_random_walk_ : (rad/s)/sqrt(s)
//   accel_bias_random_walk_: (m/s^2)/sqrt(s)
// Each is squared to convert a noise density into the variance rate Q
// expects (matches the Q*dt scaling already used in predictCovariance()).
// ---------------------------------------------------------------------------
void ESEKF::setImuNoiseParameters(float gyro_noise_density, float accel_noise_density,
                                   float gyro_bias_random_walk, float accel_bias_random_walk)
{
	// All four or none: a partially applied Q (three good blocks and one NaN)
	// is harder to diagnose than a call that visibly did nothing, and the
	// getters below let the caller read back what actually took effect.
	if (!isValidVariance(gyro_noise_density) || !isValidVariance(accel_noise_density) ||
	    !isValidVariance(gyro_bias_random_walk) || !isValidVariance(accel_bias_random_walk))
	{
		return;
	}

	gyro_noise_density_     = gyro_noise_density;
	accel_noise_density_    = accel_noise_density;
	gyro_bias_random_walk_  = gyro_bias_random_walk;
	accel_bias_random_walk_ = accel_bias_random_walk;

	float gyro_var      = gyro_noise_density_ * gyro_noise_density_;
	float accel_var     = accel_noise_density_ * accel_noise_density_;
	float gyro_bias_var = gyro_bias_random_walk_ * gyro_bias_random_walk_;
	float accel_bias_var = accel_bias_random_walk_ * accel_bias_random_walk_;

	for (int i = 0; i < 3; ++i)
	{
		Q[i][i]           = gyro_var;
		Q[3 + i][3 + i]   = accel_var;
		Q[9 + i][9 + i]   = gyro_bias_var;
		Q[12 + i][12 + i] = accel_bias_var;
	}
}

float ESEKF::getGyroNoiseDensity() const      { return gyro_noise_density_; }
float ESEKF::getAccelNoiseDensity() const     { return accel_noise_density_; }
float ESEKF::getGyroBiasRandomWalk() const    { return gyro_bias_random_walk_; }
float ESEKF::getAccelBiasRandomWalk() const   { return accel_bias_random_walk_; }

// ---------------------------------------------------------------------------
// Measurement noise getters / setters
// ---------------------------------------------------------------------------
// A measurement noise matrix is only accepted if it is finite with a strictly
// positive diagonal -- see isValidNoiseMat3(). A rejected call leaves the
// previous R in place rather than disabling the sensor silently.
void ESEKF::setAccelNoise(const float R[3][3])
{
	if (!isValidNoiseMat3(R)) return;
	std::memcpy(R_accel, R, sizeof(R_accel));
}

void ESEKF::setMagNoise(const float R[3][3])
{
	if (!isValidNoiseMat3(R)) return;
	std::memcpy(R_mag, R, sizeof(R_mag));
}

void ESEKF::setBaroNoise(float variance)
{
	// Strictly positive: kalmanUpdateScalar() requires S = h^2 P + r > 0, so a
	// zero or negative R_baro makes every barometer update fail its guard and
	// the vertical channel goes unaided with no error reported anywhere.
	if (!std::isfinite(variance) || variance <= 0.0f) return;
	R_baro = variance;
}

void ESEKF::setGPSNoise(const float R[3][3])
{
	if (!isValidNoiseMat3(R)) return;
	std::memcpy(R_gps, R, sizeof(R_gps));
}

void ESEKF::setGPSVelocityNoise(const float R[3][3])
{
	if (!isValidNoiseMat3(R)) return;
	std::memcpy(R_gps_vel, R, sizeof(R_gps_vel));
}

// EKF3-style parameterization: EK3_VELNE_M_NSE / EK3_VELD_M_NSE, given as
// standard deviations in m/s and squared into variances here.
void ESEKF::setGPSVelocityNoiseSigma(float sigma_horizontal, float sigma_vertical)
{
	if (!std::isfinite(sigma_horizontal) || !std::isfinite(sigma_vertical) ||
	    sigma_horizontal <= 0.0f || sigma_vertical <= 0.0f)
	{
		return;
	}

	std::memset(R_gps_vel, 0, sizeof(R_gps_vel));
	R_gps_vel[0][0] = sigma_horizontal * sigma_horizontal;
	R_gps_vel[1][1] = sigma_horizontal * sigma_horizontal;
	R_gps_vel[2][2] = sigma_vertical * sigma_vertical;
}

void ESEKF::getAccelNoise(float R_out[3][3]) const { std::memcpy(R_out, R_accel, sizeof(R_accel)); }
void ESEKF::getMagNoise(float R_out[3][3]) const   { std::memcpy(R_out, R_mag, sizeof(R_mag)); }
float ESEKF::getBaroNoise() const                  { return R_baro; }
void ESEKF::getGPSNoise(float R_out[3][3]) const   { std::memcpy(R_out, R_gps, sizeof(R_gps)); }
void ESEKF::getGPSVelocityNoise(float R_out[3][3]) const { std::memcpy(R_out, R_gps_vel, sizeof(R_gps_vel)); }

// A NaN threshold would make the motion gate's comparison false for every
// sample, i.e. the accelerometer would be fused as an attitude reference
// throughout every maneuver -- the exact failure the gate exists to prevent.
void ESEKF::setAccelGateThreshold(float threshold_mps2)
{
	if (!std::isfinite(threshold_mps2) || threshold_mps2 < 0.0f) return;
	accel_gate_threshold_ = threshold_mps2;
}
float ESEKF::getAccelGateThreshold() const               { return accel_gate_threshold_; }
float ESEKF::getAccelMotion() const                      { return accel_motion_; }
float ESEKF::getLastAccelNoiseInflation() const          { return accel_inflation_; }

// ---------------------------------------------------------------------------
// Reference / environment getters / setters
// ---------------------------------------------------------------------------
void ESEKF::setBaroReference(float altitude_ref)
{
	if (std::isfinite(altitude_ref)) baro_ref_ = altitude_ref;
}
float ESEKF::getBaroReference() const             { return baro_ref_; }

// <= 0 disables gating, which is a legitimate (if unwise) choice. NaN is not:
// it disables the gate by accident and reads back as "set".
void ESEKF::setInnovationGate(float nis_threshold)
{
	if (!std::isfinite(nis_threshold)) return;
	nis_gate_ = nis_threshold;
}
float ESEKF::getInnovationGate() const             { return nis_gate_; }

// Gravity is used both as the accelerometer's predicted measurement and as the
// centre of its motion gate, so a NaN here disables the gate AND poisons the
// innovation on the same update.
void ESEKF::setGravity(const Vector3f &g0)        { if (isFiniteVec(g0)) g = g0; }
Vector3f ESEKF::getGravity() const                { return g; }
void ESEKF::setMagReference(const Vector3f &mag_ref) { if (isFiniteVec(mag_ref)) mag_ref_ = mag_ref; }
void ESEKF::setMagneticDeclination(float declination_rad)
{
	if (std::isfinite(declination_rad)) mag_declination_ = declination_rad;
}
float ESEKF::getMagneticDeclination() const               { return mag_declination_; }
Vector3f ESEKF::getMagReference() const           { return mag_ref_; }
Vector3f ESEKF::getGPSReference() const           { return gps_ref_; }

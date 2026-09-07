/*
 * ESEKF.hpp
 *
 *  Created on: Aug 16, 2026
 *      Author: KAVINDU
 */

#ifndef ESEKF_ESEKF_HPP_
#define ESEKF_ESEKF_HPP_

// Maths types only -- deliberately NOT common.hpp, which drags in the STM32
// HAL. This class is pure algorithm and must stay portable; see MathTypes.hpp.
#include "MathTypes.hpp"

// ---------------------------------------------------------------------------
// BUILD AND USAGE CONTRACT -- read before integrating.
//
//  * Freestanding-friendly. No dynamic allocation, no exceptions, no RTTI, no
//    virtuals, no globals, no I/O, no <chrono>. The whole filter is one
//    self-contained object (~2.2 kB) that can live in .bss; it is safe to
//    build with -fno-exceptions -fno-rtti and to construct at file scope.
//  * Do NOT build with -ffast-math / -ffinite-math-only. Every divergence
//    guard in this filter is a std::isfinite() check, and those options tell
//    the compiler NaN and Inf cannot occur, which lets it delete the guards
//    outright. The filter would then propagate a NaN into the control loop
//    silently -- exactly the failure the guards exist to catch.
//  * Bounded work. Every loop has a compile-time bound and there is no
//    recursion, so predict()/update() are hard-real-time safe: worst case
//    equals typical. Peak stack is ~1.2 kB, in the GPS/mag update path.
//  * NOT thread-safe and not reentrant. One instance belongs to one task.
//    If sensor drivers deliver on other tasks, hand the samples over through a
//    queue rather than calling into the filter from both sides; there is no
//    internal locking and none is wanted at this rate.
//  * Single-precision throughout, except the filter clock and the geodetic
//    conversion, which are double for reasons documented at their definitions.
// ---------------------------------------------------------------------------

// 15 error states: [attitude(3), velocity(3), position(3), gyro_bias(3), accel_bias(3)]
constexpr int ESEKF_STATE_DIM = 15;

// predict() clamps dt to this. A longer gap (a scheduling stall, a dropped IMU
// batch) is not something a first-order strapdown integrator can absorb, so the
// step is truncated rather than allowed to fling the state across the map.
constexpr float ESEKF_MAX_PREDICT_DT = 0.1f; // seconds

// ---------------------------------------------------------------------------
// Aiding timeouts, in seconds of filter time (accumulated from predict()'s dt).
// Modelled on ArduPilot EKF3, which drives its nav_filter_status from the time
// since each source was last SUCCESSFULLY FUSED -- not the time since data last
// arrived. That distinction is the whole point: a receiver happily emitting
// fixes that the innovation gate rejects every time is not aiding the filter,
// and a flight controller that watches "is GPS connected" will not notice.
// ---------------------------------------------------------------------------
constexpr float ESEKF_GPS_USE_TIMEOUT   = 4.0f; // no fused GPS for this long -> using_gps goes false
constexpr float ESEKF_GPS_AID_TIMEOUT   = 5.0f; // ... this long -> position aiding lost, dead reckoning
constexpr float ESEKF_HGT_AID_TIMEOUT   = 5.0f; // no fused baro/GPS height for this long -> height aiding lost
constexpr float ESEKF_MAG_AID_TIMEOUT   = 5.0f; // no fused magnetometer for this long -> yaw unaided

// After a source has been gated out continuously for this long, the next sample
// is fused REGARDLESS of the gate. Straight from EKF3, which fuses on
// `(posCheckPassed || posTimeout || badIMUdata)`. Without this escape a single
// large transient -- a GPS jump, a hard landing, a magnetic slam -- can push the
// state so far from the measurement that every subsequent sample fails the gate
// and the source is locked out permanently. The filter then flies on dead
// reckoning while a perfectly good receiver sits there being ignored.
constexpr float ESEKF_GATE_RECOVERY_TIMEOUT = 2.0f;

// ---------------------------------------------------------------------------
// Filter status reported to the flight controller.
//
// Deliberately shaped like EKF3's nav_filter_status so an existing FC health
// check maps across. The flags a navigation mode should gate on are
// horiz_pos_valid / vert_pos_valid / horiz_vel_valid; dead_reckoning is the one
// that should trigger a failsafe.
// ---------------------------------------------------------------------------
struct ESEKFStatus
{
	bool attitude_valid;    // attitude estimate usable (filter healthy, not diverged)
	bool horiz_vel_valid;   // horizontal velocity observable from a live source
	bool vert_vel_valid;    // vertical velocity observable from a live source
	bool horiz_pos_valid;   // absolute horizontal position observable (GPS fusing)
	bool vert_pos_valid;    // vertical position observable (baro or GPS height fusing)
	bool using_gps;         // GPS was successfully FUSED recently (not merely received)
	bool dead_reckoning;    // NO absolute position source is being fused -- inertial only
	bool gps_glitching;     // the most recent GPS sample was REJECTED by the gate
	bool position_reset;    // position was snapped to a measurement recently --
	                        // the estimate is DISCONTINUOUS, see getLastPositionResetDelta()
	bool mag_aiding;        // magnetometer fusing, i.e. yaw is bounded
	bool initialized;       // initialize() completed successfully
	bool diverged;          // state or covariance went non-finite; filter is dead
};

class ESEKF {
public:
	ESEKF();

	// ---- Lifecycle ----
	void initialize(const Vector3f &accel, const Vector3f &mag, float altitude, const Vector3f &gps);
	void reset();                       // reset state and covariance to defaults
	bool isInitialized() const;

	// True once the filter has gone non-finite WITH NOTHING TO FALL BACK ON.
	// A diverged filter refuses all further predict()/update() calls; the caller
	// must reset() and re-initialize(). Poll this from the flight-control loop.
	//
	// This is deliberately hard to reach. A numerical fault after alignment is
	// REPAIRED rather than latched -- see getFaultCount() -- because latching in
	// flight leaves an aircraft with no attitude and no way back: recovery needs
	// a re-alignment, re-alignment needs a stationary airframe to read gravity
	// from, and an aircraft that has just lost its estimate is the opposite of
	// stationary. In practice this only latches on a garbage ALIGNMENT, which
	// happens on the ground, where refusing to fly is the correct answer.
	bool hasDiverged() const;

	// Number of numerical faults REPAIRED since power-on: times the state or
	// covariance went non-finite and the filter rolled back to its last good
	// nominal state and re-seeded the covariance instead of dying.
	//
	// Monotone for the life of the object -- neither reset() nor initialize()
	// clears it, because it is evidence about the airframe, not filter state.
	// It should read zero forever; anything else is a bug worth chasing, and
	// without this counter the repair would be completely silent.
	uint32_t getFaultCount() const;

	// ---- Filter status / integrity (poll this from the flight controller) ----
	ESEKFStatus getStatus() const;

	// Shorthand for getStatus().dead_reckoning. True when no absolute position
	// source has been fused for ESEKF_GPS_AID_TIMEOUT seconds, i.e. position is
	// pure inertial dead reckoning and its error is growing without bound.
	bool isDeadReckoning() const;

	// Seconds of filter time since each source was last SUCCESSFULLY fused.
	// Returns a large value if never fused. Log these.
	float getTimeSinceGPSFusion() const;
	float getTimeSinceBaroFusion() const;
	float getTimeSinceMagFusion() const;

	// Normalized innovation test ratios from the most recent attempt at each
	// source: 1.0 means the innovation sat exactly on the gate, > 1.0 means the
	// sample was rejected. EKF3 exposes the same quantity (posTestRatio etc.)
	// and it is the single most useful number to log -- a ratio creeping toward
	// 1.0 is a sensor going bad before it fails outright.
	float getGPSPosTestRatio() const;
	float getGPSVelTestRatio() const;
	float getBaroTestRatio() const;
	float getMagTestRatio() const;

	// Horizontal position uncertainty (m, 1-sigma) from P. A navigation mode
	// should refuse to engage above some threshold.
	float getHorizontalPositionError() const;
	float getVerticalPositionError() const;

	// ---- Position reset reporting (a controller MUST watch this) ----
	// When a source times out, the filter snaps position/velocity to the
	// measurement rather than fusing its way back (see resetPositionTo()). That
	// makes the estimate DISCONTINUOUS. A position or velocity controller
	// integrating this signal will see the jump as a huge instantaneous error
	// and command a violent correction unless it is told to shift its own
	// reference by the same amount. EKF3 publishes exactly this, as
	// getLastPosNorthEastReset() / getLastVelNorthEastReset().
	//
	// Returns the NED delta applied by the most recent reset, and the filter
	// time at which it happened.
	Vector3f getLastPositionResetDelta() const;
	Vector3f getLastVelocityResetDelta() const;
	float getLastPositionResetTime() const;
	// Filter time of the most recent VELOCITY reset. A velocity reset is just
	// as discontinuous as a position one -- resetVelocityTo() is called on a
	// GPS-velocity timeout -- and without its own timestamp a controller
	// reading getLastVelocityResetDelta() cannot tell a jump that happened one
	// cycle ago from one that happened at takeoff. Returns a large negative
	// sentinel if no velocity reset has occurred.
	float getLastVelocityResetTime() const;

	// Largest position jump a timeout reset is allowed to apply, in metres.
	// Modelled on EKF3's EK3_GLITCH_RAD (default 25 m). A reset larger than this
	// still happens -- refusing it would leave the filter permanently lost -- but
	// it RAISES ESEKFStatus::gps_glitching for the reporting window, on top of
	// position_reset and the delta getters, so a vehicle can distinguish a
	// routine re-centering from a teleport and failsafe instead of flying to a
	// waypoint in a frame that just moved 500 m.
	void setGlitchRadius(float radius_m);
	float getGlitchRadius() const;

	// ---- Predict / update ----
	// dt is clamped to (0, ESEKF_MAX_PREDICT_DT]; non-finite inputs are rejected.
	void predict(const Vector3f &gyro, const Vector3f &accel, float dt);

	// Each update returns false if it was rejected -- not initialized, diverged,
	// non-finite input, blocked by the motion gate, or blocked by the innovation
	// consistency gate -- and true if it was actually fused.
	bool updateAccelerometer(const Vector3f &accel);
	bool updateBarometer(float altitude);

	// Magnetometer. Takes the calibrated field as a body-frame 3-vector in
	// whatever unit the sensor reports, and fuses HEADING ONLY: one scalar,
	// with a Jacobian along the earth vertical. Roll and pitch are structurally
	// immune to it -- a magnetic disturbance can rotate the estimate about the
	// vertical and do nothing else, whatever the sensor reads.
	//
	// That is not a simplification, it is the correct model for this filter.
	// Gravity already fixes two of the three attitude axes; the magnetometer is
	// only needed for the third. Fusing the field as a 3-vector against a
	// reference fixed at alignment -- with no earth-field states to absorb the
	// difference -- feeds every local anomaly into roll and pitch as well. See
	// the note in updateMagnetometer().
	//
	// Returns false, and leaves yaw unaided, when the field is too weak or too
	// close to vertical to carry a heading.
	bool updateMagnetometer(const Vector3f &mag);

	// gps = {latitude, longitude, altitude}. UNITS: lat/lon in RADIANS, altitude
	// in metres. Most GPS drivers (NMEA, u-blox UBX-NAV-PVT) emit DEGREES --
	// use updateGPSDegrees() for those, or convert at the call site. Passing
	// degrees here produces a silent ~57x position scale error, so choose
	// deliberately. initialize()'s `gps` argument uses the same units.
	bool updateGPS(const Vector3f &gps);

	// As updateGPS(), but takes {latitude(deg), longitude(deg), altitude(m)}.
	// initialize()'s gps argument must still be given in radians.
	bool updateGPSDegrees(const Vector3f &gps_deg);

	// ---- GPS velocity (NED, m/s) ----
	// Fuse the receiver's Doppler velocity solution. Measurement model is the
	// identity on the velocity states: z = v_ned, H = [0 I 0 0 0].
	//
	// This is the most valuable measurement the filter was missing. GPS Doppler
	// velocity is far more accurate than differentiated position (typically
	// 0.05-0.1 m/s vs metres), and it is what actually makes the accelerometer
	// bias and the roll/pitch tilt observable -- the two states this filter
	// estimated worst without it. EKF3 treats GPS velocity as a primary source
	// and defaults its noise to 0.3 m/s horizontal, 0.5 m/s vertical.
	//
	// Feed this at the same rate as the position fix, from the same solution.
	bool updateGPSVelocity(const Vector3f &velocity_ned);

	// PREFERRED GPS ENTRY POINT.
	// A latitude/longitude in radians does not survive float32: one ULP at 45
	// deg latitude is 6.0e-8 rad, i.e. 0.38 m of ground distance (0.76 m at 60
	// deg). Passing lat/lon through a Vector3f therefore quantizes the position
	// measurement more coarsely than a decent GPS's own noise floor, and the
	// damage is done at the API boundary -- before the filter can do anything
	// about it. This overload keeps the angles in double all the way to the
	// (small) local NED difference, which is then exact in float.
	bool updateGPSPrecise(double latitude_rad, double longitude_rad, float altitude_m);

	// Sets the local NED origin in full precision. initialize() calls this with
	// its (float) gps argument; call it directly afterwards if you have the
	// origin in double.
	void setGPSReferencePrecise(double latitude_rad, double longitude_rad, float altitude_m);

	// ---- State resets ----
	// Hard-set a state block and re-seed its covariance, discarding the prior.
	//
	// An innovation gate can lock a healthy source out permanently: reject a
	// sample because the state disagrees, and the state drifts further, so the
	// next sample looks worse still. For a BOUNDED quantity, fusing anyway once
	// the source times out is enough to recover. For position and velocity --
	// which are unbounded and can be arbitrarily far from truth after a glitch
	// or a bad alignment -- it is not: each forced fusion only takes a small
	// bite out of the error and then restarts the timeout, so recovery takes
	// tens of seconds. Measured on a 500 m state error: 5 forced fusions in
	// 10 s, still 60 m out at the end.
	//
	// EKF3 does not fuse its way out of that. On a position timeout it calls
	// ResetPosition(): zeroStatesVarCov(7,8), snap the state to the GPS
	// measurement, and set P to the measurement noise. When the filter has
	// established that the state is wrong and the sensor is right, believe the
	// sensor. These do the same.
	void resetPositionTo(const Vector3f &position_ned, float variance);
	void resetVelocityTo(const Vector3f &velocity_ned, float variance);
	void resetVerticalPositionTo(float down_m, float variance);

	// ---- Error-state injection ----
	// Injects a 15x1 error state (dx) into the nominal state. dx is produced
	// and consumed within a single update (K * innovation -> inject -> discard),
	// so it is passed in rather than stored as a member.
	void injectErrorState(const float dx[ESEKF_STATE_DIM]);

	// ---- State getters ----
	Quaternionf getOrientation() const;
	Vector3f getEulerAngles() const;      // roll, pitch, yaw (rad) -> returned as x,y,z
	Vector3f getVelocity() const;
	Vector3f getPosition() const;         // local NED position, meters
	Vector3f getGPSPosition() const;      // {latitude(rad), longitude(rad), altitude(m)} -- inverse of updateGPS()'s conversion, using getGPSReference() as origin
	Vector3f getGyroBias() const;
	Vector3f getAccelBias() const;
	Vector3f getAngularRate() const;      // body-frame angular rate (rad/s), bias-corrected, from the most recent predict()
	float getAltitude() const;            // current altitude estimate (m), derived from position, not the raw last barometer reading

	// ---- State setters ----
	void setOrientation(const Quaternionf &q0);
	void setVelocity(const Vector3f &v0);
	void setPosition(const Vector3f &p0);
	void setGyroBias(const Vector3f &b0);
	void setAccelBias(const Vector3f &b0);

	// ---- Covariance getters / setters ----
	void getCovariance(float P_out[ESEKF_STATE_DIM][ESEKF_STATE_DIM]) const;
	void setCovariance(const float P_in[ESEKF_STATE_DIM][ESEKF_STATE_DIM]);
	float getStateVariance(int index) const; // P[index][index], index in [0, 15)

	// ---- Process noise getters / setters ----
	void setProcessNoiseGyro(float sigma2);
	void setProcessNoiseAccel(float sigma2);
	void setProcessNoiseGyroBias(float sigma2);
	void setProcessNoiseAccelBias(float sigma2);
	void getProcessNoise(float Q_out[ESEKF_STATE_DIM][ESEKF_STATE_DIM]) const;

	// Preferred way to set Q: from IMU datasheet noise-density figures rather
	// than hand-picked sigma^2 values. Squares each input internally and
	// writes the result onto Q's diagonal.
	//   gyro_noise_density   : (rad/s)/sqrt(Hz)   -- "angular random walk"
	//   accel_noise_density  : (m/s^2)/sqrt(Hz)   -- "velocity random walk"
	//   gyro_bias_random_walk: (rad/s)/sqrt(s)    -- "bias instability"
	//   accel_bias_random_walk: (m/s^2)/sqrt(s)
	void setImuNoiseParameters(float gyro_noise_density, float accel_noise_density,
	                            float gyro_bias_random_walk, float accel_bias_random_walk);
	float getGyroNoiseDensity() const;
	float getAccelNoiseDensity() const;
	float getGyroBiasRandomWalk() const;
	float getAccelBiasRandomWalk() const;

	// ---- Measurement noise getters / setters ----
	void setAccelNoise(const float R[3][3]);

	// Magnetometer noise, as a covariance in the sensor's own field units --
	// the same units updateMagnetometer() is fed. Only the diagonal average is
	// used: the fusion is scalar (heading), so the field-unit variance is
	// converted to a heading variance by dividing by the squared horizontal
	// field strength, which is the correct conversion and also makes a weak
	// horizontal field automatically trusted less. Keeping the setter in field
	// units means a caller tunes against a number the datasheet quotes.
	void setMagNoise(const float R[3][3]);
	void setBaroNoise(float variance);
	void setGPSNoise(const float R[3][3]);
	void setGPSVelocityNoise(const float R[3][3]);
	// Convenience: isotropic horizontal + separate vertical, as EKF3 parameterizes
	// it (EK3_VELNE_M_NSE / EK3_VELD_M_NSE). Arguments are standard deviations
	// in m/s; squared internally.
	void setGPSVelocityNoiseSigma(float sigma_horizontal, float sigma_vertical);
	void getAccelNoise(float R_out[3][3]) const;
	void getMagNoise(float R_out[3][3]) const;
	float getBaroNoise() const;
	void getGPSNoise(float R_out[3][3]) const;
	void getGPSVelocityNoise(float R_out[3][3]) const;

	// Motion gate for updateAccelerometer(): the update is skipped whenever
	// |accel| deviates from local gravity magnitude by more than this many
	// m/s^2 (the vehicle is assumed to be under real acceleration, not just
	// gravity, so the reading is not a trustworthy "up" reference). Default
	// is set in the constructor; tune tighter to reject more motion, looser
	// to keep correcting through vibration/maneuvering. Set <= 0 to disable
	// the gate, which also disables the noise inflation below.
	void setAccelGateThreshold(float threshold_mps2);
	float getAccelGateThreshold() const;

	// | |accel| - |g| | from the most recent updateAccelerometer() call,
	// whether it was fused or gated out. This is the quantity the gate above is
	// compared against, so it is the number to tune the gate FROM: hover the
	// vehicle, read this, and set the gate just above the value it settles at.
	float getAccelMotion() const;

	// Variance multiplier R_accel was scaled by for the most recently FUSED
	// accelerometer sample; sqrt() it for the sigma multiplier. Always >= 1.
	//
	// A bare threshold is a cliff: a sample just inside it is trusted in full
	// and one just outside is thrown away entirely, so a vehicle sitting near
	// the gate flips between full aiding and none from cycle to cycle. Instead
	// the measurement noise is inflated smoothly with how far into the gate the
	// sample sits, so a marginal sample is de-weighted rather than trusted or
	// discarded outright, and the gate only draws the line where the reading
	// stops being a gravity reference at all.
	float getLastAccelNoiseInflation() const;

	// ---- Reference / environment getters / setters ----
	void setGravity(const Vector3f &g0);
	Vector3f getGravity() const;
	void setMagReference(const Vector3f &mag_ref);
	Vector3f getMagReference() const;

	// ---- Magnetic declination: this is what makes yaw ABSOLUTE ----
	// Angle from TRUE north to MAGNETIC north, radians, positive east.
	//
	// Without it this filter's yaw is only relative: initialize() previously
	// derived mag_ref_ by rotating the first magnetometer sample through the
	// attitude it had just computed from that same sample, so the reference was
	// self-consistent by construction and carried no absolute information. The
	// magnetometer then held yaw wherever initialization happened to put it, and
	// any initial heading error was permanent -- an autonomous vehicle would fly
	// straight, confident, and in a rotated frame.
	//
	// With a declination set, initialize() instead builds the NED reference
	// field as {H*cos(dec), H*sin(dec), D}, so its horizontal component points
	// at true magnetic north. This is what EKF3's alignMagStateDeclination()
	// does. Yaw then means true heading.
	//
	// Get the value from a World Magnetic Model lookup for your operating
	// location (ArduPilot ships AP_Declination; it varies from roughly -25 to
	// +25 degrees across the populated world, and is about -2 deg in Sri Lanka).
	// MUST be set before initialize() to take effect on the initial attitude.
	void setMagneticDeclination(float declination_rad);
	float getMagneticDeclination() const;
	Vector3f getGPSReference() const; // raw GPS reading at inilialize() (origin of the local NED frame)

	// ---- Barometric altitude datum ----
	// The barometer is fused as a RELATIVE altitude sensor: updateBarometer()
	// forms its residual against (altitude - baro_ref_), the climb above the
	// initialization point, which is directly comparable to -position_NED_.z
	// whatever datum the barometer itself uses (pressure altitude referenced to
	// 1013.25 hPa, QNH, MSL...).
	//
	// This datum must NOT be conflated with the GPS altitude datum in
	// gps_ref_.z: barometric pressure altitude and GPS ellipsoidal height
	// routinely differ by tens to hundreds of metres, and feeding that offset
	// into the residual injects it straight into the vertical state.
	// Set automatically by initialize() from its `altitude` argument.
	void setBaroReference(float altitude_ref);
	float getBaroReference() const;

	// ---- Innovation consistency gating ----
	// An update is rejected when its normalized innovation squared
	// (NIS = nu^T S^-1 nu) exceeds this threshold, keeping a single bad sample
	// -- GPS multipath, a magnetic transient, a barometric pressure spike --
	// out of the state. Units: chi-square with dim(z) degrees of freedom; the
	// default is a deliberately loose gate. Set <= 0 to disable gating.
	void setInnovationGate(float nis_threshold);
	float getInnovationGate() const;

private:
	// State vector
	Quaternionf q;              // Orientation quaternion (nominal state)
	Vector3f velocity_NED_;     // Velocity, NED frame
	Vector3f position_NED_;     // Position, local NED frame (origin at gps_ref_)
	Vector3f gyro_bias;         // Gyroscope bias
	Vector3f accel_bias;        // Accelerometer bias

	// Error covariance matrix
	float P[ESEKF_STATE_DIM][ESEKF_STATE_DIM]; // 15x15 error covariance matrix

	// Process noise covariance
	float Q[ESEKF_STATE_DIM][ESEKF_STATE_DIM]; // 15x15 process noise covariance matrix

	// Measurement noise covariance
	float R_accel[3][3]; // 3x3 accelerometer measurement noise covariance
	float R_mag[3][3];   // 3x3 magnetometer measurement noise covariance
	float R_baro;        // Barometer measurement noise variance
	float R_gps[3][3];     // 3x3 GPS position measurement noise covariance
	float R_gps_vel[3][3]; // 3x3 GPS velocity (NED) measurement noise covariance
	float accel_gate_threshold_; // m/s^2, see setAccelGateThreshold()
	float accel_motion_;         // | |accel| - |g| | from the last update attempt
	float accel_inflation_;      // R_accel multiplier used by the last fused sample

	// IMU continuous-time noise parameters (datasheet / Allan-variance figures).
	// Kept alongside Q so the values used to derive it are inspectable/reusable,
	// e.g. if Q is rebuilt after changing dt characteristics.
	float gyro_noise_density_;       // (rad/s)/sqrt(Hz)
	float accel_noise_density_;      // (m/s^2)/sqrt(Hz)
	float gyro_bias_random_walk_;    // (rad/s)/sqrt(s)
	float accel_bias_random_walk_;   // (m/s^2)/sqrt(s)

	// Reference vectors
	Vector3f g;           // Gravity vector (NED)
	Vector3f angular_rate_; // last bias-corrected body-frame angular rate (rad/s), cached by predict()
	Vector3f mag_ref_;    // Reference magnetic field vector (NED)
	// Full-precision copy of the local NED origin's angular coordinates. gps_ref_
	// below keeps the float32 view for the existing getGPSReference() API, but
	// every conversion uses these.
	double gps_ref_lat_rad_;
	double gps_ref_lon_rad_;

	Vector3f gps_ref_;    // Raw GPS reading at inilialize() — origin of the local NED frame.
	                       // updateGPS() must convert incoming raw GPS against this before
	                       // treating it as a NED position (e.g. via an ECEF/LLA -> local NED
	                       // conversion), it is not itself in NED.

	// Bookkeeping
	bool initialized_;
	bool diverged_;     // latched only when a fault has no last-good state to fall back on

	// ---- Numerical fault recovery ----
	// The last NOMINAL state that passed checkFinite(), and nothing else: 16
	// floats, refreshed on every healthy predict()/update(). Snapshotting P too
	// would be 900 more floats copied at the loop rate for no benefit -- a
	// covariance that has gone bad is not worth restoring, and resetting it to
	// the initial diagonal is both cheaper and the honest statement of what the
	// filter knows after a fault (see checkFinite()).
	Quaternionf last_good_q_;
	Vector3f    last_good_velocity_;
	Vector3f    last_good_position_;
	Vector3f    last_good_gyro_bias_;
	Vector3f    last_good_accel_bias_;
	bool        have_last_good_;  // false until the first healthy check after alignment
	uint32_t    fault_count_;     // repairs performed; see getFaultCount()
	float dt_last_;
	float altitude_;    // last raw barometer altitude reading (not the filter's estimate -- see getAltitude())
	float baro_ref_;    // barometer datum -- the raw reading at initialize(). See setBaroReference().
	float nis_gate_;    // innovation gate threshold. See setInnovationGate().
	float mag_declination_; // radians, positive east. See setMagneticDeclination().
	float glitch_radius_;   // m, see setGlitchRadius()

	// Reset reporting
	Vector3f last_pos_reset_delta_;
	Vector3f last_vel_reset_delta_;
	double   last_pos_reset_t_;
	double   last_vel_reset_t_;
	bool     last_pos_reset_exceeded_glitch_; // |horizontal reset| > glitch_radius_

	// ---- Aiding bookkeeping ----
	// filter_time_ accumulates dt from predict(), giving a monotonic clock with
	// no dependency on the caller having a wall clock. Everything below is
	// stamped against it. Sentinel -1000 means "never fused".
	//
	// DOUBLE, deliberately, and this is not premature caution. In float32 the
	// clock stops advancing once it outgrows its own increment: at 200 Hz
	// (dt = 5 ms) the ULP of a float exceeds dt at 2^24 * 5 ms, and rounding
	// starts biting long before that -- past ~36 h of continuous running,
	// filter_time_ += dt is a no-op. Every aiding timeout below is a difference
	// against this clock, so a frozen clock means "time since last fusion"
	// freezes too: the filter would report GPS aiding forever, dead_reckoning
	// never, and no failsafe would ever fire. That is a silent, unrecoverable
	// integrity failure on exactly the systems that stay powered longest.
	// One double add per predict() (~100 cycles on a soft-float M4F, 0.002% of
	// a 200 Hz budget) buys a clock good for longer than any vehicle's life.
	double filter_time_;
	double last_gps_pos_fuse_t_;
	double last_gps_vel_fuse_t_;
	double last_baro_fuse_t_;
	double last_mag_fuse_t_;
	double last_gps_pos_reject_t_; // last time a GPS position sample was GATED OUT
	double last_gps_vel_reject_t_;
	double last_baro_reject_t_;
	double last_mag_reject_t_;

	// Most recent normalized test ratios (innovation^2 / gate). 1.0 = on the gate.
	float pos_test_ratio_;
	float vel_test_ratio_;
	float hgt_test_ratio_;
	float mag_test_ratio_;

	// Set by kalmanUpdate()/kalmanUpdateScalar() so the caller-facing update
	// wrappers can record pass/fail and the test ratio without duplicating the
	// gate arithmetic.
	float last_test_ratio_;
	bool  last_update_gated_;

	// ---- Helper functions: quaternion operations (stateless -- const) ----
	Quaternionf quaternionMultiply(const Quaternionf &q1, const Quaternionf &q2) const;
	Quaternionf quaternionConjugate(const Quaternionf &q) const;
	Quaternionf quaternionNormalize(const Quaternionf &q) const;
	Quaternionf quaternionFromEuler(float roll, float pitch, float yaw) const;
	// Exact exponential map SO(3) <- R^3. Replaces the previous first-order
	// [1, dtheta/2] form, which truncates the rotation angle by |dtheta|^3/12
	// per step and cannot represent a large error-state correction.
	Quaternionf quaternionFromRotationVector(const Vector3f &delta_theta) const;
	Vector3f eulerFromQuaternion(const Quaternionf &q) const;
	Vector3f rotateVector(const Quaternionf &q, const Vector3f &v) const;
	void quaternionToRotationMatrix(const Quaternionf &q, float R_out[3][3]) const;

	// ---- Helper functions: linear algebra (stateless -- const) ----
	Vector3f vectorAdd(const Vector3f &a, const Vector3f &b) const;
	Vector3f vectorSub(const Vector3f &a, const Vector3f &b) const;
	Vector3f vectorScale(const Vector3f &a, float s) const;
	Vector3f crossProduct(const Vector3f &a, const Vector3f &b) const;
	float vectorNorm(const Vector3f &a) const;
	void skewSymmetric(const Vector3f &v, float M_out[3][3]) const;

	void matrixMultiply3x3(const float A[3][3], const float B[3][3], float C_out[3][3]) const;
	bool matrixInverse3x3(const float A[3][3], float A_inv_out[3][3]) const;

	// ---- Filter internals ----
	// P <- F P F^T + Q*dt, with F the discrete error-state transition matrix
	// I + A*dt linearized at the MID-INTERVAL attitude R_mid -- the same point
	// the nominal velocity integration uses, so nominal and error-state
	// propagation stay consistent. omega and f arrive already bias-corrected.
	//
	// F is never MATERIALIZED. It is I plus five sparse blocks, so the product
	// is applied as in-place row and column operations on P. Building F densely
	// and multiplying costs 900 bytes of stack for F, another 900 for the F*P
	// temporary and ~6750 multiply-accumulates; the structured form costs a
	// handful of scalars and ~1000. On an STM32 running FreeRTOS with 1-2 kB
	// task stacks that 1.8 kB is the difference between running and a stack
	// overflow, and predict() is the highest-rate call in the filter.
	void predictCovariance(const Vector3f &omega, const Vector3f &f,
	                        const float R_mid[3][3], float dt);

	// Vector (3-element) and scalar measurement updates. Both use the Joseph
	// form and both return false if the innovation gate rejected the sample.
	// force_fuse bypasses the innovation gate. Used only by the timeout recovery
	// path -- see ESEKF_GATE_RECOVERY_TIMEOUT.
	bool kalmanUpdate(const float innovation[3], const float H[3][ESEKF_STATE_DIM], const float R[3][3],
	                   bool force_fuse = false);
	bool kalmanUpdateScalar(float innovation, int state_index, float h, float r,
	                         bool force_fuse = false);

	// Scalar update whose Jacobian is non-zero on exactly one 3-state block,
	// H = [0.. h3 ..0] starting at `first_state`. The magnetometer's heading
	// measurement is this shape: one number, observing a direction in the
	// attitude block. kalmanUpdateScalar() cannot express it (its H has a single
	// non-zero ELEMENT) and the 3-vector kalmanUpdate() would need a 3x3 inverse
	// to fuse one scalar. Joseph form, O(N^2), 60 bytes of scratch.
	bool kalmanUpdateScalarVec3(float innovation, int first_state, const float h3[3],
	                             float r, bool force_fuse = false);

	// True if `last_reject_t` shows this source has been continuously gated out
	// for longer than ESEKF_GATE_RECOVERY_TIMEOUT, meaning the next sample must
	// be fused regardless of the gate to break the lockout.
	bool gateRecoveryDue(double last_fuse_t, double last_reject_t) const;

	// True once a source has gone longer than `timeout` without a successful
	// fusion, i.e. it is no longer aiding the filter at all.
	bool sourceTimedOut(double last_fuse_t, float timeout) const;

	// Zeros rows and columns [first, last] of P and sets their diagonal to
	// `variance` -- the covariance half of a state reset. Equivalent to EKF3's
	// zeroStatesVarCov() followed by seeding the diagonal.
	void resetStateBlockCovariance(int first, int last, float variance);

	void resetCovarianceDefaults(); // sets P to the default diagonal (used by both the constructor and inilialize())

	// Folds the error state into the nominal state AND applies the ESEKF
	// covariance reset P <- G P G^T. Every update path goes through this rather
	// than calling injectErrorState() directly.
	void injectAndResetCovariance(const float dx[ESEKF_STATE_DIM]);

	void symmetrizeCovariance();  // P <- 0.5 * (P + P^T)

	// Health check run after every propagation and every fusion. Returns true
	// while the filter is healthy, and on the way through captures the nominal
	// state as the rollback point. On a fault it REPAIRS -- rolls the nominal
	// state back and re-seeds P -- and only latches diverged_ when there is no
	// rollback point yet, which means the fault came from alignment itself.
	bool checkFinite();
};

#endif /* ESEKF_ESEKF_HPP_ */

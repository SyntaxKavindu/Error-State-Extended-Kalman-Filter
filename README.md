# Error-State Extended Kalman Filter (ESEKF)

A 15-state error-state extended Kalman filter for inertial navigation on
microcontrollers. It fuses a 6-axis IMU with a magnetometer, barometer and GPS
(position **and** Doppler velocity) to produce attitude, velocity, position and
sensor-bias estimates — together with the integrity flags a flight controller
needs in order to decide whether those estimates can be trusted.

Written for an STM32F7 flight controller, but the filter itself carries no
vendor, HAL or RTOS dependency: it is freestanding C++11 that builds and runs
anywhere, including a host-side test harness.

```
IMU (gyro + accel) ──► predict()  ──┐
accelerometer ──────► update…()  ──┤
magnetometer ───────► update…()  ──┼──► ESEKF ──► attitude / velocity / position
barometer ──────────► update…()  ──┤              gyro & accel bias
GPS pos + velocity ─► update…()  ──┘              ESEKFStatus (integrity flags)
```

---

## Table of contents

- [Why error-state](#why-error-state)
- [Files](#files)
- [Conventions: frames, units, sign](#conventions-frames-units-sign)
- [Quick start](#quick-start)
- [State vector](#state-vector)
- [Measurement models](#measurement-models)
- [Integrity reporting](#integrity-reporting)
- [Position and velocity resets](#position-and-velocity-resets)
- [Robustness](#robustness)
- [Tuning](#tuning)
- [API reference](#api-reference)
- [Building](#building)
- [Constraints and non-goals](#constraints-and-non-goals)
- [References](#references)

---

## Why error-state

A direct EKF estimates the state itself. That is awkward for attitude: a
quaternion has four parameters for three degrees of freedom, so its covariance
is structurally singular and the normalization constraint has to be re-imposed
after every update.

An **error-state** filter splits the problem in two:

- the **nominal state** — quaternion, velocity, position, biases — is integrated
  by the full non-linear strapdown mechanization, with no linearization at all;
- the **error state** — a 15-element vector of small deviations from the nominal
  — is what the Kalman filter actually carries a covariance for.

Attitude error is parameterized as a 3-element rotation vector, so the
covariance is minimal and non-singular. The error state stays close to zero,
which keeps the linearization valid, and it is injected into the nominal state
and reset to zero at the end of every fusion — so it never has a chance to grow
large enough for the first-order approximation to break down.

This filter follows that structure throughout, and follows ArduPilot's EKF3
closely on the operational details — timeouts, gating, resets, status flags —
because those are where an estimator that is mathematically correct still gets a
vehicle lost.

## Files

| File | Contents |
|---|---|
| `ESEKF.hpp` / `ESEKF.cpp` | The filter. Everything else is support. |
| `MathTypes.hpp` | Umbrella header — pulls in the three types below plus `<cmath>` / `<cstring>`. |
| `Vector3f.hpp` | 3-vector: arithmetic, `length()`, `dot()`. |
| `Quaternionf.hpp` | Quaternion: Hamilton product, `rotate()`, `conjugate()`. |
| `Matrix3f.hpp` | 3×3 matrix: products, `transposed()`, `scaled()`, `identity()`. |

The maths types are deliberately kept in their own header rather than in a
project-wide `common.hpp`. Including a board header to get `Vector3f` is what
makes an algorithm unbuildable anywhere but the board it was written on.

## Conventions: frames, units, sign

Getting these wrong is the most common integration failure, and most of them
fail *silently* — a plausible-looking estimate in the wrong frame.

| Quantity | Convention |
|---|---|
| Navigation frame | **NED** — x North, y East, z **Down** |
| Body frame | **FRD** — x Forward, y Right, z Down |
| Gravity | `g_ned = {0, 0, +9.80665}` — down is positive |
| Attitude | Hamilton quaternion, body→NED, scalar-first `(w, x, y, z)` |
| Euler angles | 3-2-1 (roll, pitch, yaw), radians, returned in a `Vector3f` as `(x=roll, y=pitch, z=yaw)` |
| Angular rate | rad/s, body frame |
| Acceleration | m/s², body frame, **specific force** (a level, stationary vehicle reads `{0, 0, −9.81}`) |
| Velocity / position | m/s and m, NED, origin at the `initialize()` GPS fix |
| Altitude | metres, **up-positive** at the API; internally `−position.z` |
| Latitude / longitude | **radians** in `updateGPS()` and `initialize()`; use `updateGPSDegrees()` for degrees |
| Magnetic field | any consistent unit — only the direction is used |
| Declination | radians, positive **east** |

**On GPS units.** `updateGPS()` takes radians. Nearly every GPS driver — NMEA,
u-blox `UBX-NAV-PVT` — emits degrees, and passing those to the radian API is a
silent ~57× position scale error. Three entry points exist so the choice is
deliberate:

```cpp
filter.updateGPS(Vector3f(lat_rad, lon_rad, alt_m));      // radians
filter.updateGPSDegrees(Vector3f(lat_deg, lon_deg, alt_m)); // degrees
filter.updateGPSPrecise(lat_rad, lon_rad, alt_m);         // preferred — double
```

`updateGPSPrecise()` is the one to use. A latitude in radians does not survive
`float32`: one ULP at 45° latitude is 6.0e-8 rad, which is 0.38 m of ground
distance (0.76 m at 60°). Passing lat/lon through a `Vector3f` therefore
quantizes the measurement more coarsely than a decent receiver's own noise
floor, and the damage is done at the API boundary before the filter can do
anything about it. The `double` overload keeps the angles in full precision all
the way to the (small) local NED difference, which is then exact in `float`.

## Quick start

```cpp
#include "ESEKF.hpp"

ESEKF filter;   // safe at file scope: trivially destructible, no allocation

void setup()
{
    // 1. Tune BEFORE initializing — declination in particular must be set
    //    first for it to affect the initial heading.
    filter.setMagneticDeclination(-0.035f);   // rad, positive east (~-2 deg)

    filter.setImuNoiseParameters(
        1.0e-3f,   // gyro noise density,        (rad/s)/sqrt(Hz)
        2.0e-2f,   // accel noise density,       (m/s^2)/sqrt(Hz)
        3.0e-4f,   // gyro bias random walk,     (rad/s)/sqrt(s)
        3.0e-3f);  // accel bias random walk,    (m/s^2)/sqrt(s)

    filter.setGPSVelocityNoiseSigma(0.3f, 0.5f);  // m/s, horizontal / vertical

    // 2. Align, with the vehicle STATIONARY. Roll and pitch come from
    //    gravity, heading from the tilt-compensated magnetic field.
    //    initialize() can fail — check it.
    filter.initialize(accel_sample,      // m/s^2, body FRD
                      mag_sample,        // body FRD, any unit
                      baro_altitude_m,   // m, up-positive
                      Vector3f(lat_rad, lon_rad, gps_alt_m));

    if (!filter.isInitialized()) {
        // Bad alignment sample — retry with good data. Do not fly.
    }
}

// IMU rate, e.g. 200 Hz
void onImuSample(const Vector3f &gyro, const Vector3f &accel, float dt)
{
    filter.predict(gyro, accel, dt);
}

// Slower, whenever each sensor delivers
void onAccel(const Vector3f &a) { filter.updateAccelerometer(a); }
void onMag  (const Vector3f &m) { filter.updateMagnetometer(m); }
void onBaro (float alt_m)       { filter.updateBarometer(alt_m); }

void onGps(double lat_rad, double lon_rad, float alt_m, const Vector3f &vel_ned)
{
    filter.updateGPSPrecise(lat_rad, lon_rad, alt_m);
    filter.updateGPSVelocity(vel_ned);   // fuse both, from the same solution
}

// Control loop
void control()
{
    const ESEKFStatus st = filter.getStatus();

    if (st.diverged) { /* reset() + initialize(); no estimate available */ }
    if (!st.attitude_valid) { /* attitude cannot be trusted */ }
    if (st.dead_reckoning)  { /* no absolute position source — failsafe */ }
    if (st.position_reset)  { /* shift the controller's reference, see below */ }

    const Vector3f euler = filter.getEulerAngles();  // roll, pitch, yaw (rad)
    const Vector3f vel   = filter.getVelocity();     // NED, m/s
    const Vector3f pos   = filter.getPosition();     // local NED, m
    const Vector3f rate  = filter.getAngularRate();  // body, rad/s, de-biased
}
```

Every `update*()` returns `bool`: `true` if the sample was actually fused,
`false` if it was rejected — not initialized, diverged, non-finite input, gated
out by the motion gate, or gated out by the innovation consistency check. The
return value is worth logging; a source that is being received but never fused
is not aiding the filter at all.

### Alignment requirements

`initialize()` takes its sample entirely on trust — there is no prior to gate it
against — so whatever arrives becomes the attitude, the NED origin and the
barometric datum. A bad sample here does not produce an estimate that converges
away; it produces a filter that is confidently wrong for the whole flight.

- The vehicle **must be stationary**. Roll and pitch are derived from the
  assumption that the accelerometer reads gravity alone.
- The accelerometer magnitude must be within ±50% of `|g|`, or alignment is
  refused. This catches a dead sensor (reading zero looks exactly like a level
  vehicle) or a saturated one, both of which would otherwise align just as
  confidently as a good part.
- A zero magnetic field is *not* an error: yaw simply starts at zero and the
  local frame is not north-aligned until some other source constrains it. Roll
  and pitch remain fully observable without a magnetometer.
- `setMagneticDeclination()` must be called **before** `initialize()` to affect
  the initial heading.
- On failure, `isInitialized()` stays `false` and every `predict()`/`update*()`
  is a no-op. **Check it.**

## State vector

The error state is 15 elements:

| Index | Error state | Nominal state | Model |
|---|---|---|---|
| 0–2 | attitude error `δθ` (rad) | quaternion `q` | rotation vector, body-frame (right) multiplication |
| 3–5 | velocity error `δv` (m/s) | `velocity_NED_` | strapdown integration of specific force |
| 6–8 | position error `δp` (m) | `position_NED_` | trapezoidal integration of velocity |
| 9–11 | gyro bias error `δb_g` (rad/s) | `gyro_bias` | random walk |
| 12–14 | accel bias error `δb_a` (m/s²) | `accel_bias` | random walk |

### Propagation

`predict()` runs the nominal state through a strapdown mechanization and the
covariance through `P ← F P Fᵀ + Q·dt`:

- **Attitude** advances by the **exact exponential map** `q ← q ⊗ Exp(ω·dt)`.
  The common first-order `[1, δθ/2]` form truncates the rotation angle by
  `|δθ|³/12` per step and cannot represent a large error-state correction.
- **Velocity** is integrated with the **mid-interval attitude**. Rotating
  specific force by the attitude at either endpoint is first-order accurate but
  carries half a step of misalignment, worth `R(ω × f)·dt²/2` of velocity error
  per step — a systematic, coherent error under sustained rotation. The midpoint
  costs one extra quaternion product and makes the integration second-order.
- **Position** uses the trapezoidal average of the old and new velocity.
- **`F` is never materialized.** It is the identity plus five sparse blocks, so
  the product is applied as in-place row and column operations on `P`: about
  1000 multiply-accumulates and a few dozen bytes of stack, against ~6750 MACs
  and 1.8 kB for the dense form. On a task with a 1–2 kB stack that difference
  is the difference between running and overflowing. The structured result is
  numerically identical to the dense product to float round-off.
- `dt` is clamped to `ESEKF_MAX_PREDICT_DT` (0.1 s). A longer gap — a scheduling
  stall, a dropped IMU batch — is not something a first-order integrator can
  absorb, so the step is truncated rather than allowed to fling the state across
  the map.

Note on `Q`: the noise-input matrix `G` is taken as the identity, which is
correct only because every `Q` block is isotropic (`σ²·I`). Accelerometer noise
enters the error dynamics as `−R·n_a`, whose contribution is
`R(σ²I)Rᵀ = σ²I` — rotation-invariant. Making any `Q` block anisotropic requires
reinstating `G` and rotating it.

## Measurement models

All updates use the **Joseph form**, which stays symmetric and positive
semi-definite under round-off where the short `(I − KH)P` form does not.

| Source | `z` | Jacobian | Notes |
|---|---|---|---|
| Accelerometer | `Rᵀ(−g) + b_a` | `[[v]ₓ, 0, 0, 0, I]` | tilt reference; motion-gated |
| Magnetometer | heading (scalar) | `h₃` on the attitude block only | **heading only** |
| Barometer | `−p_z` (scalar) | `H[8] = −1` | relative to the alignment datum |
| GPS position | `p_ned` | `[0, 0, I, 0, 0]` | WGS84 local-tangent conversion |
| GPS velocity | `v_ned` | `[0, I, 0, 0, 0]` | Doppler solution — exact, no linearization |

### The magnetometer fuses heading only, and this matters

The textbook construction fuses the field as a 3-vector against a reference
fixed at alignment. That is wrong for a filter shaped like this one, because
there are no earth-field states to absorb the difference between the reference
and the field actually present. Every discrepancy — a local anomaly, a steel
gate, current through a nearby harness, flying somewhere with a different
inclination than where the vehicle was switched on — is then fused as a full
3-axis attitude error, and two thirds of that correction lands in **roll and
pitch**.

A vehicle told it is banked when it is level does exactly the wrong thing, and
because the accelerometer fights the error back, the symptom is not a clean
failure: it is attitude that leans toward magnetic clutter and recovers slowly.
Very hard to attribute in flight.

Geomagnetism carries no information about tilt anyway — gravity already fixes
two axes, and the magnetometer is only needed for the third. So this filter
fuses exactly that: one scalar, the heading discrepancy, with a Jacobian that
can only rotate the estimate about the earth vertical. **Roll and pitch are
structurally immune to the magnetometer, whatever it reads.** This is what
ArduPilot's `fuseEulerYaw()` does.

The measurement noise is set in the sensor's own field units
(`setMagNoise()` — the units the datasheet quotes) and converted internally to a
heading variance by dividing by the squared horizontal field strength, which
also makes a weak horizontal field automatically trusted less. Near-vertical or
vanishing fields are refused outright: `updateMagnetometer()` returns `false`
and yaw is left unaided, which the status flags then report honestly.

### Declination is what makes yaw absolute

Without a declination, the natural reference construction —
`mag_ref_ = rotate(q, mag)` — rotates the first sample through the attitude just
derived *from that same sample*. The result is self-consistent by construction:
the first innovation is identically zero, and stays zero for any initial yaw
whatsoever. The magnetometer then holds yaw wherever initialization happened to
put it and supplies no absolute information — an initial heading error becomes
permanent and unobservable, and an autonomous vehicle flies straight, confident,
and in a rotated frame.

With a declination set, `initialize()` builds the NED reference field as
`{H·cos(dec), H·sin(dec), D}`, so its horizontal component points at true
magnetic north and **yaw means true heading**. Get the value from a World
Magnetic Model lookup for your operating area (ArduPilot ships
`AP_Declination`); it ranges roughly ±25° across the populated world.

### The accelerometer gate is graded, not binary

`updateAccelerometer()` uses the accelerometer as a gravity reference, which is
only valid when the vehicle is not under real acceleration. The sample is
skipped when `| |a| − |g| |` exceeds `setAccelGateThreshold()` (default 1.0 m/s²).

Inside the gate, trust is graded rather than binary: `R` is inflated by
`1 + 9·(motion/gate)²`, so a sample sitting on the boundary is charged about
3.2× its nominal sigma instead of being taken at full confidence one cycle and
discarded the next. A bare threshold makes a vehicle hovering near the limit
chatter between full aiding and none, which shows up as attitude that twitches
with throttle.

To tune the gate, hover the vehicle, read `getAccelMotion()` — reported whether
the sample was fused or not, so it still works while everything is being
rejected — and set the threshold just above where it settles.

### GPS velocity is the most valuable measurement here

Doppler velocity is accurate to ~0.05–0.1 m/s where differentiated position is
accurate to metres, and because the velocity error dynamics couple to attitude
through `−R[f]ₓ` and to accelerometer bias through `−R`, it is what actually
makes **tilt and accelerometer bias observable**. Those are the two states this
filter estimates worst without it. Feed `updateGPSVelocity()` at the same rate
as the position fix, from the same solution.

### Geodetic conversion

`updateGPSPrecise()` converts to the local tangent plane using the WGS84
ellipsoid with the **meridional** radius of curvature for north and the
**prime-vertical** radius (times `cos(lat)`) for east. Using the equatorial
radius for both — a common shortcut — puts a ~0.7% scale error on the north
axis. Longitude differences are wrapped into `[−π, π]`, so a track crossing the
antimeridian does not produce a 2π jump in the east coordinate.

## Integrity reporting

This is the interface an autonomous vehicle actually flies on. Position,
velocity and attitude are only meaningful read *together with* these flags,
because an inertial filter keeps producing smooth, plausible, confidently-wrong
output long after its absolute references have gone away. That is not a defect
to be fixed — it is what inertial navigation *does* — and the only defence is to
say so out loud.

```cpp
struct ESEKFStatus {
    bool attitude_valid;   // attitude usable (healthy, not diverged)
    bool horiz_vel_valid;  // horizontal velocity observable from a live source
    bool vert_vel_valid;   // vertical velocity observable
    bool horiz_pos_valid;  // absolute horizontal position observable
    bool vert_pos_valid;   // vertical position observable (baro or GPS height)
    bool using_gps;        // GPS successfully FUSED recently
    bool dead_reckoning;   // NO absolute position source — inertial only
    bool gps_glitching;    // the most recent GPS sample was REJECTED
    bool position_reset;   // position snapped — the estimate is DISCONTINUOUS
    bool mag_aiding;       // magnetometer fusing, i.e. yaw is bounded
    bool initialized;      // initialize() completed successfully
    bool diverged;         // state/covariance non-finite; filter is dead
};
```

Gate navigation modes on `horiz_pos_valid` / `vert_pos_valid` /
`horiz_vel_valid`. `dead_reckoning` is the one that should trigger a failsafe.

**Every flag is driven by the time since a source was last *successfully
fused*** — never by whether data is arriving. A receiver happily emitting fixes
that the innovation gate rejects on every sample is contributing exactly
nothing, and a health check that watches "GPS connected" or satellite count will
report that everything is fine.

The timing comes from `filter_time_`, a monotonic clock accumulated from
`predict()`'s `dt`, so the filter needs no wall clock from the caller. It is a
`double` deliberately: in `float32` at 200 Hz the clock stops advancing once it
outgrows its own increment — past roughly 36 hours of continuous running,
`filter_time_ += dt` becomes a no-op. Every timeout is a difference against that
clock, so a frozen clock means the filter reports GPS aiding forever,
`dead_reckoning` never, and no failsafe ever fires. One double add per predict
(~100 cycles on a soft-float M4F, 0.002% of a 200 Hz budget) buys a clock good
for longer than any vehicle's service life.

| Constant | Value | Meaning |
|---|---|---|
| `ESEKF_GPS_USE_TIMEOUT` | 4.0 s | no fused GPS for this long → `using_gps` false |
| `ESEKF_GPS_AID_TIMEOUT` | 5.0 s | → position aiding lost, dead reckoning |
| `ESEKF_HGT_AID_TIMEOUT` | 5.0 s | → height aiding lost |
| `ESEKF_MAG_AID_TIMEOUT` | 5.0 s | → yaw unaided |
| `ESEKF_GATE_RECOVERY_TIMEOUT` | 2.0 s | gated out this long → force-fuse the next sample |

### Test ratios

`getGPSPosTestRatio()`, `getGPSVelTestRatio()`, `getBaroTestRatio()` and
`getMagTestRatio()` return the normalized innovation test ratio from the most
recent attempt at each source: `1.0` means the innovation sat exactly on the
gate, `> 1.0` means rejected. **Log these.** A ratio creeping toward 1.0 is a
sensor going bad before it fails outright — the single most useful diagnostic
the filter exposes.

`getHorizontalPositionError()` and `getVerticalPositionError()` give 1-sigma
uncertainties in metres from `P`; a navigation mode should refuse to engage
above some threshold.

## Position and velocity resets

An innovation gate can lock a healthy source out permanently: reject a sample
because the state disagrees, and the state drifts further, so the next sample
looks worse still. Two mechanisms break that loop.

**Forced fusion.** After a source has been gated out continuously for
`ESEKF_GATE_RECOVERY_TIMEOUT` (2 s), the next sample is fused regardless of the
gate. This is enough for a *bounded* quantity.

**State reset.** For position and velocity it is not: they are unbounded and can
be arbitrarily far from truth after a glitch or a bad alignment, and each forced
fusion only takes a small bite out of the error before restarting the timeout.
Measured on a 500 m state error: five forced fusions in 10 s, still 60 m out at
the end. So on a position or height timeout the filter does what EKF3 does —
snap the state to the measurement and re-seed its covariance from the
measurement noise. When the filter has established that the state is wrong and
the sensor is right, believe the sensor.

**This makes the estimate discontinuous, and a controller must be told.** A
position or velocity controller integrating this signal sees the jump as a huge
instantaneous error and commands a violent correction unless it shifts its own
reference by the same amount:

```cpp
if (filter.getStatus().position_reset) {
    const Vector3f d = filter.getLastPositionResetDelta();  // NED metres
    setpoint_ned = setpoint_ned + d;        // move the reference with the frame
}
```

`getLastVelocityResetDelta()` is the velocity equivalent, and
`getLastPositionResetTime()` / `getLastVelocityResetTime()` timestamp each
against the filter clock, so a controller can tell a jump one cycle old from one
at takeoff. Both flags are reported for a full aiding window so a slow control
loop cannot miss them.

**Glitch radius.** `setGlitchRadius()` (default 25 m, EKF3's `EK3_GLITCH_RAD`)
separates a routine few-metre re-centering after a brief outage from a 500 m
teleport. A reset larger than the radius still happens — refusing it would leave
the filter permanently lost — but it additionally raises `gps_glitching`, so a
vehicle can failsafe instead of flying to a waypoint in a frame that just moved.

## Robustness

**Divergence is repaired, not latched.** Every `predict()` and every fusion runs
a health check. On the way through it snapshots the nominal state — 16 floats,
not `P` — as a rollback point. If the state or covariance goes non-finite, the
filter rolls back to that snapshot, re-seeds `P` to its default diagonal, and
increments `getFaultCount()`.

Latching in flight would leave an aircraft with no attitude and no way back:
recovery needs a re-alignment, re-alignment needs a stationary airframe to read
gravity from, and an aircraft that has just lost its estimate is the opposite of
stationary. `hasDiverged()` therefore only latches when there is *no* rollback
point — which means the fault came from alignment itself, on the ground, where
refusing to fly is the correct answer.

`getFaultCount()` is monotone for the life of the object; neither `reset()` nor
`initialize()` clears it, because it is evidence about the airframe rather than
filter state. **It should read zero forever** — anything else is a bug worth
chasing, and without the counter the repair would be completely silent.

**Negative variances are handled by magnitude, not sign.** A diagonal entry of
−1e−30 on a state whose variance is 1e−9 is float32 round-off in the Joseph
product, not divergence, and latching on it would kill a healthy filter for a
rounding error. The test is relative to the size of the covariance: round-off is
clamped to zero (as EKF3 and PX4's ekf2 do in `ConstrainVariances()`) and
re-inflated by `Q` on the next predict; anything meaningfully negative still
latches.

**Non-finite inputs are rejected everywhere.** `dt <= 0.0f` alone does not reject
a NaN — every comparison against NaN is false in IEEE-754, so a NaN would walk
straight through such a guard and turn the whole state into NaN with no recovery
path. Finiteness is tested explicitly on every input.

**Every setter validates.** These values arrive from a parameter store, a
ground-station link or an auto-tune routine, and a single NaN written into `R`,
`Q` or gravity propagates into `P` on the next update. A rejected setter leaves
the previous, known-good tuning in place — always the safer failure. Noise
matrices are additionally required to have a strictly positive diagonal, since a
zero or negative variance makes `S = HPHᵀ + R` singular or indefinite.
`setImuNoiseParameters()` applies all four values or none: a partially applied
`Q` is harder to diagnose than a call that visibly did nothing.

## Tuning

### Defaults

| Parameter | Default | Setter |
|---|---|---|
| Gyro noise density | 1.0e−3 (rad/s)/√Hz | `setImuNoiseParameters()` |
| Accel noise density | 2.0e−2 (m/s²)/√Hz | `setImuNoiseParameters()` |
| Gyro bias random walk | 3.0e−4 (rad/s)/√s | `setImuNoiseParameters()` |
| Accel bias random walk | 3.0e−3 (m/s²)/√s | `setImuNoiseParameters()` |
| Accelerometer noise | 0.05 (diagonal) | `setAccelNoise()` |
| Magnetometer noise | 0.05 (diagonal, field units²) | `setMagNoise()` |
| Barometer noise | 1.0 m² | `setBaroNoise()` |
| GPS position noise | 1.0 m² (diagonal) | `setGPSNoise()` |
| GPS velocity noise | σ = 0.3 / 0.5 m/s (h / v) | `setGPSVelocityNoiseSigma()` |
| Accel motion gate | 1.0 m/s² | `setAccelGateThreshold()` |
| Innovation gate (NIS) | 25.0 | `setInnovationGate()` |
| Glitch radius | 25 m | `setGlitchRadius()` |
| Magnetic declination | 0 rad | `setMagneticDeclination()` |

Initial covariance diagonal: attitude 0.05 rad², velocity 0.25 (m/s)², position
1.0 m², gyro bias 1e−4 (rad/s)², accel bias 1e−3 (m/s²)².

The IMU and GPS-velocity defaults follow a low-cost MEMS part and EKF3's
`EK3_VELNE_M_NSE` / `EK3_VELD_M_NSE` respectively. They are starting points, not
recommendations.

### Where to start

1. **`setImuNoiseParameters()` first.** Take the four figures from your IMU's
   datasheet or an Allan-variance run rather than hand-picking `σ²` values —
   they are squared internally onto `Q`'s diagonal. This is the single highest-
   leverage tuning step.
2. **Measurement noise from datasheets** where possible. GPS position and
   velocity accuracy are usually quoted directly; barometer noise is easiest to
   measure by logging a stationary unit for a minute.
3. **`setAccelGateThreshold()` from `getAccelMotion()`**, as described above.
4. **The innovation gate last.** The default 25 is deliberately permissive
   (χ² at 99.9% is 16.3 for 3 DOF, 10.8 for 1 DOF), so it only catches grossly
   inconsistent samples and never healthy transients. Tighten it only with test
   ratio logs in hand. Set `≤ 0` to disable gating entirely.

Tuning values are *not* filter state: `initialize()` and `reset()` leave every
noise parameter, gate, gravity vector and declination alone, so a re-alignment
does not discard them. The *references* are a different matter — `reset()`
clears the barometric datum and `initialize()` re-derives both it and the NED
origin from the new alignment sample, which is the point of re-aligning.

## API reference

### Lifecycle

```cpp
void     initialize(const Vector3f &accel, const Vector3f &mag,
                    float altitude, const Vector3f &gps);
void     reset();                  // state and covariance back to defaults
bool     isInitialized() const;
bool     hasDiverged() const;
uint32_t getFaultCount() const;
```

### Predict and update

```cpp
void predict(const Vector3f &gyro, const Vector3f &accel, float dt);

bool updateAccelerometer(const Vector3f &accel);
bool updateMagnetometer (const Vector3f &mag);
bool updateBarometer    (float altitude);
bool updateGPS          (const Vector3f &gps);         // lat/lon RADIANS
bool updateGPSDegrees   (const Vector3f &gps_deg);     // lat/lon degrees
bool updateGPSPrecise   (double lat_rad, double lon_rad, float alt_m);  // preferred
bool updateGPSVelocity  (const Vector3f &velocity_ned);
```

### State access

```cpp
Quaternionf getOrientation() const;
Vector3f    getEulerAngles() const;   // roll, pitch, yaw (rad) as x, y, z
Vector3f    getVelocity()    const;   // NED, m/s
Vector3f    getPosition()    const;   // local NED, m
Vector3f    getGPSPosition() const;   // {lat(rad), lon(rad), alt(m)}
Vector3f    getGyroBias()    const;
Vector3f    getAccelBias()   const;
Vector3f    getAngularRate() const;   // body, rad/s, de-biased, from last predict
float       getAltitude()    const;   // m, from the estimate — not the raw baro

void setOrientation(const Quaternionf &q0);
void setVelocity   (const Vector3f &v0);
void setPosition   (const Vector3f &p0);
void setGyroBias   (const Vector3f &b0);
void setAccelBias  (const Vector3f &b0);
```

`getAngularRate()` is the bias-corrected body rate cached by the last
`predict()` — this is what a rate controller should read, rather than the raw
gyro.

### Status and integrity

```cpp
ESEKFStatus getStatus() const;
bool  isDeadReckoning() const;

float getTimeSinceGPSFusion()  const;   // seconds of filter time
float getTimeSinceBaroFusion() const;
float getTimeSinceMagFusion()  const;

float getGPSPosTestRatio() const;       // 1.0 = on the gate, >1.0 = rejected
float getGPSVelTestRatio() const;
float getBaroTestRatio()   const;
float getMagTestRatio()    const;

float getHorizontalPositionError() const;   // m, 1-sigma
float getVerticalPositionError()   const;

Vector3f getLastPositionResetDelta() const; // NED metres
Vector3f getLastVelocityResetDelta() const;
float    getLastPositionResetTime()  const;
float    getLastVelocityResetTime()  const;

float getAccelMotion()              const; // | |a| - |g| |, last attempt
float getLastAccelNoiseInflation()  const; // R multiplier, last fused sample
```

### Resets

```cpp
void resetPositionTo        (const Vector3f &position_ned, float variance);
void resetVelocityTo        (const Vector3f &velocity_ned, float variance);
void resetVerticalPositionTo(float down_m, float variance);
void injectErrorState       (const float dx[ESEKF_STATE_DIM]);
```

### Covariance, process and measurement noise

```cpp
void  getCovariance(float P_out[15][15]) const;
void  setCovariance(const float P_in[15][15]);
float getStateVariance(int index) const;        // index in [0, 15)

void setProcessNoiseGyro(float sigma2);
void setProcessNoiseAccel(float sigma2);
void setProcessNoiseGyroBias(float sigma2);
void setProcessNoiseAccelBias(float sigma2);
void getProcessNoise(float Q_out[15][15]) const;
void setImuNoiseParameters(float gyro_nd, float accel_nd,
                           float gyro_brw, float accel_brw);

void setAccelNoise(const float R[3][3]);
void setMagNoise  (const float R[3][3]);        // sensor field units
void setBaroNoise (float variance);
void setGPSNoise  (const float R[3][3]);
void setGPSVelocityNoise(const float R[3][3]);
void setGPSVelocityNoiseSigma(float sigma_h, float sigma_v);  // m/s
```

Matching getters exist for all of these.

### Reference and environment

```cpp
void     setGravity(const Vector3f &g0);
Vector3f getGravity() const;
void     setMagReference(const Vector3f &mag_ref);
Vector3f getMagReference() const;
void     setMagneticDeclination(float declination_rad);   // before initialize()
float    getMagneticDeclination() const;
Vector3f getGPSReference() const;
void     setGPSReferencePrecise(double lat_rad, double lon_rad, float alt_m);
void     setBaroReference(float altitude_ref);
float    getBaroReference() const;

void  setAccelGateThreshold(float threshold_mps2);
void  setInnovationGate(float nis_threshold);   // <= 0 disables gating
void  setGlitchRadius(float radius_m);
```

The barometric datum must **not** be conflated with the GPS altitude datum:
pressure altitude and GPS ellipsoidal height routinely differ by tens to
hundreds of metres, and feeding that offset into the residual injects it
straight into the vertical state. `initialize()` sets both correctly and
separately.

## Building

Drop the five files into your project and compile `ESEKF.cpp`. There is nothing
to configure and nothing to link.

```sh
# Host build / test harness
g++ -std=c++11 -O2 -Wall -Wextra -c ESEKF.cpp

# Cortex-M7F with hardware single-precision FPU
arm-none-eabi-g++ -std=c++11 -O2 \
    -mcpu=cortex-m7 -mfpu=fpv5-sp-d16 -mfloat-abi=hard \
    -fno-exceptions -fno-rtti -c ESEKF.cpp
```

**Requirements**

- A **C++11 freestanding** implementation with `<cmath>` and `<cstring>`.
- IEEE-754 binary32 `float`. Enforced by `static_assert`.

**Do not build with `-ffast-math` or `-ffinite-math-only`.** Every divergence
guard in this filter is an `std::isfinite()` check, and those options tell the
compiler that NaN and Inf cannot occur — which licenses it to delete the guards
outright. The filter would then feed a NaN straight into the control loop with
`hasDiverged()` still returning `false`. `ESEKF.cpp` contains an `#error` guard
on `__FAST_MATH__`, so this fails the build rather than the flight.

**Footprint and timing**

| Property | Value |
|---|---|
| `sizeof(ESEKF)` | 2352 bytes — fits in `.bss`, safe at file scope |
| Peak stack | ~1.2 kB, in the GPS/mag update path |
| Dynamic allocation | none |
| Exceptions / RTTI / virtuals / globals / I/O | none |
| Worst-case execution | bounded — every loop has a compile-time bound, no recursion |

Because there is no recursion and no data-dependent iteration, `predict()` and
`update*()` are hard-real-time safe: worst case equals typical.

`static_assert`s enforce that `ESEKF` stays trivially destructible (a static
instance needs no `atexit` slot) and trivially copyable (a snapshot can be
`memcpy`'d out for logging).

## Constraints and non-goals

- **Not thread-safe, not reentrant.** One instance belongs to one task. If
  sensor drivers deliver on other tasks, hand samples over through a queue
  rather than calling in from both sides; there is no internal locking and none
  is wanted at this rate.
- **Single precision throughout**, except the filter clock and the geodetic
  conversion, which are `double` for the reasons documented above.
- **No earth-field states.** The magnetometer is heading-only by design (see
  above); estimating the earth field the way EKF3 does would need six more
  states and is out of scope.
- **No airspeed, optical-flow, range-finder or beacon fusion**, and no multi-IMU
  or multi-GPS voting. This is one filter over one sensor set.
- **No wind or terrain estimation.**
- **Calibration is the caller's job.** The filter expects a calibrated
  magnetometer (hard- and soft-iron corrected) and accelerometer; it estimates
  *biases* as states, but does not perform calibration.

## References

The operational behaviour follows ArduPilot's EKF3 closely, and the derivations
follow the standard error-state references:

- ArduPilot **EKF3** (`AP_NavEKF3`) — `nav_filter_status`, aiding timeouts, gate
  recovery, `ResetPosition()`, `EK3_GLITCH_RAD`, `fuseEulerYaw()`,
  `alignMagStateDeclination()`, `ConstrainVariances()`.
- PX4 **ecl/ekf2** — variance constraint handling.
- J. Solà, *Quaternion kinematics for the error-state Kalman filter* (2017).
- P. D. Groves, *Principles of GNSS, Inertial, and Multisensor Integrated
  Navigation Systems*, 2nd ed.

---

`ESEKF.hpp` and `ESEKF.cpp` carry substantially more detail than this README —
every non-obvious choice is documented at the point where it was made, usually
with the failure mode that motivated it. Read them before changing anything.

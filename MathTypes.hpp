/*
 * MathTypes.hpp
 *
 *  Created on: Sep 6, 2026
 *      Author: KAVINDU
 *
 * The vector/quaternion/matrix types and the freestanding C++ library pieces
 * that go with them -- and NOTHING ELSE. No HAL, no board headers, no
 * peripherals.
 *
 * This exists so the parts of this tree that are pure algorithm -- ESEKF,
 * AccelerometerCalibrator, CompassCalibrator, LevelCalibrator -- can be lifted
 * into another project, another MCU, or a host-side test harness without
 * carrying an STM32F7 dependency they never use. Including common.hpp for
 * Vector3f alone pulled stm32f7xx_hal.h in behind it, which made every one of
 * those classes unbuildable anywhere else.
 *
 * common.hpp is now this header plus the HAL, so the drivers -- which genuinely
 * do need HAL_GPIO_WritePin() and friends -- are unaffected.
 *
 * Requirements on the toolchain: a C++11 freestanding implementation with
 * <cmath> and <cstring>. No dynamic allocation, no exceptions, no RTTI.
 */

#ifndef COMMON_MATHTYPES_HPP_
#define COMMON_MATHTYPES_HPP_

#include <stdint.h>
#include <stddef.h>
#include <cmath>   // sqrtf, fabsf, atan2f, cosf, sinf, acosf
#include <cstring> // memset, memcpy

#include "Vector3f.hpp"
#include "Quaternionf.hpp"
#include "Matrix3f.hpp"

#endif /* COMMON_MATHTYPES_HPP_ */

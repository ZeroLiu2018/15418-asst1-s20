#include <stdio.h>
#include <algorithm>
#include <math.h>
#include "CMU418intrin.h"
#include "logger.h"
using namespace std;

void absSerial(float *values, float *output, int N) {
  for (int i = 0; i < N; i++) {
    float x = values[i];
    if (x < 0) {
      output[i] = -x;
    } else {
      output[i] = x;
    }
  }
}

// implementation of absolute value using 15418 instrinsics
void absVector(float *values, float *output, int N) {
  __cmu418_vec_float x;
  __cmu418_vec_float result;
  __cmu418_vec_float zero = _cmu418_vset_float(0.f);
  __cmu418_mask maskAll, maskIsNegative, maskIsNotNegative;

  //  Note: Take a careful look at this loop indexing.  This example
  //  code is not guaranteed to work when (N % VECTOR_WIDTH) != 0.
  //  Why is that the case?
  for (int i = 0; i < N; i += VECTOR_WIDTH) {

    // All ones
    maskAll = _cmu418_init_ones();

    // All zeros
    maskIsNegative = _cmu418_init_ones(0);

    // Load vector of values from contiguous memory addresses
    _cmu418_vload_float(x, values + i, maskAll);               // x = values[i];

    // Set mask according to predicate
    _cmu418_vlt_float(maskIsNegative, x, zero, maskAll);     // if (x < 0) {

    // Execute instruction using mask ("if" clause)
    _cmu418_vsub_float(result, zero, x, maskIsNegative);      //   output[i] = -x;

    // Inverse maskIsNegative to generate "else" mask
    maskIsNotNegative = _cmu418_mask_not(maskIsNegative);     // } else {

    // Execute instruction ("else" clause)
    _cmu418_vload_float(result, values + i, maskIsNotNegative); //   output[i] = x; }

    // Write results back to memory
    _cmu418_vstore_float(output + i, result, maskAll);
  }
}

// Accepts an array of values and an array of exponents
// For each element, compute values[i]^exponents[i] and clamp value to
// 4.18.  Store result in outputs.
// Uses iterative squaring, so that total iterations is proportional
// to the log_2 of the exponent
void clampedExpSerial(float *values, int *exponents, float *output, int N) {
  for (int i = 0; i < N; i++) {
    float x = values[i];
    float result = 1.f;
    int y = exponents[i];
    float xpower = x;
    while (y > 0) {
      if (y & 0x1) {
        result *= xpower;
        if (result > 4.18f) {
          result = 4.18f;
          break;
        }
      }
      xpower = xpower * xpower;
      y >>= 1;
    }
    output[i] = result;
  }
}

void clampedExpVector(float *values, int *exponents, float *output, int N) {
  // TODO: Implement your vectorized version of clampedExpSerial here
  __cmu418_vec_int y{};
  __cmu418_vec_float xpower{};
  auto maskAll = _cmu418_init_ones();
  auto zero = _cmu418_vset_int(0);
  auto one = _cmu418_vset_int(1);
  auto four_18 = _cmu418_vset_float(4.18f);
  int aligN = N - N % VECTOR_WIDTH;
  for (int i = 0; i < aligN; i += VECTOR_WIDTH) {
    auto result = _cmu418_vset_float(1.0f);
    _cmu418_vload_float(xpower, values + i, maskAll);
    _cmu418_vload_int(y, exponents + i, maskAll);
    auto maskCur = _cmu418_init_ones();
    while (_cmu418_cntbits(maskCur)) {
      __cmu418_vec_int oddMaskInt{};
      _cmu418_vbitand_int(oddMaskInt, y, one, maskCur);
      auto oddMask = _cmu418_init_ones(0);
      _cmu418_vgt_int(oddMask, oddMaskInt, zero, maskCur);
      _cmu418_vmult_float(result, result, xpower, oddMask);
      auto larger4_18 = _cmu418_init_ones(0);
      _cmu418_vgt_float(larger4_18, result, four_18, oddMask);
      _cmu418_vset_float(result, 4.18f, larger4_18);

      auto not_larger4_18 = _cmu418_mask_not(larger4_18);
      maskCur = _cmu418_mask_and(maskCur, not_larger4_18);

      _cmu418_vmult_float(xpower, xpower, xpower, maskCur);
      _cmu418_vshiftright_int(y, y, one, maskCur);
      _cmu418_vneq_int(maskCur, y, zero, maskCur);
    }
    _cmu418_vstore_float(output + i, result, maskAll);
  }
  clampedExpSerial(values + aligN, exponents + aligN, output + aligN, N - aligN);
}

float arraySumSerial(float *values, int N) {
  float sum = 0;
  for (int i = 0; i < N; i++) {
    sum += values[i];
  }

  return sum;
}

// Assume N % VECTOR_WIDTH == 0
// Assume VECTOR_WIDTH is a power of 2
float arraySumVector(float *values, int N) {
  // TODO: Implement your vectorized version here
  auto sumVec = _cmu418_vset_float(0.0f);
  __cmu418_vec_float v{};
  auto mask_all = _cmu418_init_ones();
  int aligN = N - N % VECTOR_WIDTH;
  for (int i = 0; i < aligN; i += VECTOR_WIDTH) {
    _cmu418_vload_float(v, values + i, mask_all);
    _cmu418_vadd_float(sumVec, sumVec, v, mask_all);
  }
  float sum[VECTOR_WIDTH];
  _cmu418_vstore_float(sum, sumVec, mask_all);
  return arraySumSerial(sum, 8) + arraySumSerial(values + aligN, N - aligN);
}

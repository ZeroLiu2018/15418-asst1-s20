#include <algorithm>

// Generate random data
void initRandom(float *values, int N) {
  for (int i = 0; i < N; i++) {
    // random input values
    values[i] = 0.f + 1.99f * static_cast<float>(rand()) / RAND_MAX;
  }
}

// Generate data that gives high relative speedup
void initGood(float *values, int N) {
  float interval = 1.0 / N;
  for (int i = 0; i < N; i += 8) {
    auto v = 0.f + interval * i;
    for (int j = 0; j < 8; ++j) {
      values[i + 1] = v;
    }
  }
}

// Generate data that gives low relative speedup
void initBad(float *values, int N) {
  for (int i = 0; i < N; ++i) {
    values[i] = i % 8 == 0 ? 0.f : 1.f;
  }
}


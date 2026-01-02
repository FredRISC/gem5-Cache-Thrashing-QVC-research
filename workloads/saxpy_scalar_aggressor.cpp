/**
 * project_1b_aggressor_scalar.cpp
 *
 * Scalar Aggressor (Core 1) - v2 (Looping)
 *
 * Runs the *identical* SAXPY kernel in a loop, but as
 * simple scalar code. This is the "control" experiment.
 * It will have the same memory footprint but a "smoother"
 * memory access pattern.
 */
#include <iostream>
#include <vector>
#include <cstdint>

// 4MB of floats (4 * 1024 * 1024 bytes)
// 1M floats = 4MB
#define VECTOR_SIZE (1024 * 1024)
#define NUM_LOOPS 100 // Loop to create continuous contention

// Use float for SAXPY
std::vector<float> x(VECTOR_SIZE);
std::vector<float> y(VECTOR_SIZE);
float a = 3.14159f;

/**
 * Runs SAXPY using standard scalar C++ code.
 * This creates "smooth" memory traffic.
 */
void saxpy_scalar() {
    for (size_t i = 0; i < VECTOR_SIZE; ++i) {
        y[i] = a * x[i] + y[i];
    }
}

int main() {
    // Initialize arrays (to prevent CoW optimizations/CoW page fault)
    for (size_t i = 0; i < VECTOR_SIZE; ++i) {
        x[i] = static_cast<float>(i);
        y[i] = 1.0f;
    }

    std::cout << "Aggressor (Core 1): Starting SCALAR stream  (infinite loop)." << std::endl;

    // This is the main contention loop
    while(1) {
        saxpy_scalar();
    }
    
    std::cout << "Aggressor (Core 1): SCALAR stream complete." << std::endl;

    return 0;
}


/**
 * project_1b_aggressor_vector.cpp
 *
 * Vector Aggressor (Core 1) - v3 (Infinite Loop)
 *
 * Runs a SAXPY kernel in an infinite loop to create continuous,
 * "bursty" L2 cache pollution.
 *
 * This process is now a "slave" and will be terminated
 * by the victim (Core 0) calling m5_exit().
 */
#include <iostream>
#include <vector>
#include <cstdint>
#include <riscv_vector.h> // For RISC-V vector intrinsics

// 8MB of floats (4 * 1024 * 2048 bytes)
// 2Mi floats = 8MiB
#define VECTOR_SIZE (1024 * 2048)

// Use float for SAXPY
std::vector<float> x(VECTOR_SIZE);
std::vector<float> y(VECTOR_SIZE);
float a = 3.14159f;

/**
 * Runs SAXPY using RISC-V Vector Intrinsics.
 * This creates "bursty" memory traffic.
 */
void saxpy_vector() {
    size_t n = VECTOR_SIZE;
    size_t vl; // Vector Length
    
    for (size_t i = 0; i < n; i += vl) {
        // Set vector length for this iteration
        vl = __riscv_vsetvl_e32m8(n - i);

        // Load vector x (vle.v)
        vfloat32m8_t vx = __riscv_vle32_v_f32m8(&x[i], vl);
        
        // Load vector y (vle.v)
        vfloat32m8_t vy = __riscv_vle32_v_f32m8(&y[i], vl);

        // Compute y = a*x + y (vfmacc.vf)
        vy = __riscv_vfmacc_vf_f32m8(vy, a, vx, vl);

        // Store vector y (vse.v)
        __riscv_vse32_v_f32m8(&y[i], vy, vl);
    }
}

int main() {
    // Initialize arrays (to prevent CoW optimizations/CoW page fault)
    for (size_t i = 0; i < VECTOR_SIZE; ++i) {
        x[i] = static_cast<float>(i);
        y[i] = 1.0f;
    }

    std::cout << "Aggressor (Core 1): Starting VECTOR stream (infinite loop)." << std::endl;
    
    // This is the main contention loop
    // It will run forever until Core 0 calls m5_exit().
    while (true) {
        saxpy_vector();
    }

    // This line is never reached
    std::cout << "Aggressor (Core 1): VECTOR stream complete." << std::endl;
    
    return 0;
}


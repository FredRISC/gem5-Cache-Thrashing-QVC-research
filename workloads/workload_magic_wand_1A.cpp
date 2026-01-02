/**
 * workload_magic_wand.cpp (CORRECTED)
 *
 * This C++ program is the complete workload for Project 1A, Experiment A1.
 * It models a realistic, interactive "Magic Wand" tool workflow in an image editor
 * to demonstrate temporal cache pollution on a single-core RISC-V processor.
 *
 * The experiment is controlled by a compile-time macro `RUN_INTERFERENCE_PHASE`.
 */

// Compile with a RISC-V C++ Toolchain: riscv64-unknown-elf-g++ -O3 -march=rv64gcv -static -DGEM5_M5OPS_H -I. -o workload_magic_wand_cpp workload_magic_wand.cpp


#include <iostream>
#include <vector>
#include <cstdlib>
#include <string>
#include <cmath> // For std::abs

// RISC-V Vector intrinsics are C-compatible and can be included in C++
#include <riscv_vector.h>

// Include the gem5 m5ops header. This allows the program to control the simulator.
#ifdef GEM5_M5OPS_H
#include "gem5/m5ops.h"
#endif

// --- Configuration Constants based on the Proposal ---

// The full image is 2048x2048 pixels. Each pixel is 1 byte (grayscale).
// This creates a 4MB dataset for the VECTOR AGGRESSOR.
constexpr int IMAGE_WIDTH = 2048;
constexpr int IMAGE_HEIGHT = 2048;
constexpr long FULL_IMAGE_SIZE = IMAGE_WIDTH * IMAGE_HEIGHT;

// Define a smaller tile for the scalar task.
// Total working set = 128KB (64KB image tile + 64KB visited map), which fits in the 256KB L2.
constexpr int TILE_WIDTH = 256;
constexpr int TILE_HEIGHT = 256;

// The metadata for the scalar task's frontier queue.
constexpr int FRONTIER_QUEUE_SIZE = 1024;


// A meaningful computation to simulate work done per-pixel in the scalar task.
void update_selection_stats(long pixel_count, double& mean, double& m2, unsigned char new_pixel_value) {
    double delta = new_pixel_value - mean;
    mean += delta / pixel_count;
    double delta2 = new_pixel_value - mean;
    m2 += delta * delta2;
}


// --- Data Structures for the "Magic Wand" ---

struct Coordinate {
    int x;
    int y;
};

// --- Scalar "Victim" Task  ---
void magic_wand_select(std::vector<unsigned char>& full_image, std::vector<bool>& visited_map_tile, 
                       std::vector<Coordinate>& queue, Coordinate tile_offset, Coordinate start_pos_in_tile, 
                       unsigned char target_color, int tolerance) {
    std::cout << "Phase: Starting SCALAR work (Magic Wand Select on Tile)." << std::endl;

    long pixel_count = 0;
    double mean_color = 0.0;
    double m2_color = 0.0;

    // Reset metadata for a fresh selection
    std::fill(visited_map_tile.begin(), visited_map_tile.end(), false);
    queue.clear();
    
    int queue_head = 0;

    // Start with the initial pixel (coordinates are now local to the tile)
    queue.push_back(start_pos_in_tile);
    visited_map_tile[start_pos_in_tile.y * TILE_WIDTH + start_pos_in_tile.x] = true;

    while (queue_head < queue.size()) {
        Coordinate current_in_tile = queue[queue_head++];
        pixel_count++;

        Coordinate neighbors_in_tile[4] = {
            {current_in_tile.x, current_in_tile.y - 1}, {current_in_tile.x, current_in_tile.y + 1},
            {current_in_tile.x - 1, current_in_tile.y}, {current_in_tile.x + 1, current_in_tile.y}
        };

        for (int i = 0; i < 4; ++i) {
            Coordinate n_tile = neighbors_in_tile[i];
            
            // Boundary check is now against TILE dimensions
            if (n_tile.x >= 0 && n_tile.x < TILE_WIDTH && n_tile.y >= 0 && n_tile.y < TILE_HEIGHT) {
                long tile_index = n_tile.y * TILE_WIDTH + n_tile.x;
                
                // Convert tile coordinates to global image coordinates for pixel value lookup
                long global_index = (tile_offset.y + n_tile.y) * IMAGE_WIDTH + (tile_offset.x + n_tile.x);

                if (!visited_map_tile[tile_index] && std::abs(full_image[global_index] - target_color) <= tolerance) {
                    update_selection_stats(pixel_count + queue.size() - queue_head, mean_color, m2_color, full_image[global_index]);
                    visited_map_tile[tile_index] = true;
                    if (queue.size() < FRONTIER_QUEUE_SIZE) {
                        queue.push_back(n_tile);
                    }
                }
            }
        }
    }
    std::cout << "Phase: Finished SCALAR work. Selected " << queue.size() << " pixels within the tile." << std::endl;
}


// --- Vector "Aggressor" Task (Unchanged) ---
void adjust_brightness(std::vector<unsigned char>& image_data, int brightness) {
    std::cout << "Phase: Starting VECTOR work (Adjust Brightness)." << std::endl;
    size_t n = image_data.size();
    unsigned char* data_ptr = image_data.data();

    for (size_t i = 0; i < n; ) {
        size_t vl = __riscv_vsetvl_e8m8(n - i);
        vuint8m8_t pixels = __riscv_vle8_v_u8m8(&data_ptr[i], vl);
        vuint8m8_t adjusted_pixels = __riscv_vsaddu_vx_u8m8(pixels, brightness, vl);
        __riscv_vse8_v_u8m8(&data_ptr[i], adjusted_pixels, vl);
        i += vl;
    }
    std::cout << "Phase: Finished VECTOR work." << std::endl;
}


// --- Main Orchestrator ---
int main() {
    // --- EXPERIMENT CONTROL ---
    //#define RUN_INTERFERENCE_PHASE

    // --- Memory Allocation ---
    std::vector<unsigned char> full_image(FULL_IMAGE_SIZE);
    
    // The visited map is now sized for the TILE
    std::vector<bool> visited_map_tile(TILE_WIDTH * TILE_HEIGHT);
    std::vector<Coordinate> frontier_queue;
    frontier_queue.reserve(FRONTIER_QUEUE_SIZE);

    // Create a predictable image: a dark background with a large, bright square in the middle
    std::fill(full_image.begin(), full_image.end(), 10);
    for (int y = 500; y < 1500; ++y) {
        for (int x = 500; x < 1500; ++x) {
            full_image[y * IMAGE_WIDTH + x] = 100; // Bright square
        }
    }
    
    // Define the tile location and the start pixel within it
    Coordinate tile_top_left = {900, 900}; // Place the tile's corner inside the bright square
    Coordinate start_pixel_in_tile = {0, 0}; // Start flood-fill from the tile's bottom left

    // ===================================================================
    // THE EXPERIMENT
    // ===================================================================

    // 1. Prime Phase (Scalar on Tile)
    magic_wand_select(full_image, visited_map_tile, frontier_queue, tile_top_left, start_pixel_in_tile, 100, 5);

    #ifdef RUN_INTERFERENCE_PHASE
    // 2. Pollution Phase (Vector on Full Image)
    adjust_brightness(full_image, 50);
    #endif

    // 3. Measurement Phase (Scalar on Tile)
    #ifdef GEM5_M5OPS_H
    m5_reset_stats(0, 0);
    #endif

    // The target color is changed to 150 because adjust_brightness added 50 to the original 100.
    magic_wand_select(full_image, visited_map_tile, frontier_queue, tile_top_left, start_pixel_in_tile, 100, 5);

    #ifdef GEM5_M5OPS_H
    m5_dump_stats(0, 0);
    #endif

    return 0;
}
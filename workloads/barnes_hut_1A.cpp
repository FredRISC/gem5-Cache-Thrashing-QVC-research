/**
 * workload_barnes_hut.cpp (Final Version)
 *
 * This C++ program implements a Barnes-Hut N-body simulation for Project 1A, Exp A2.
 * The scalar task (force calculation) is a realistic, cache-sensitive workload.
 * The vector task (position update) is a streaming, cache-polluting workload.
 */

/* Compile: Run ./workload.sh or Compile for specific modes as per below commands

# INTERFERENCE MODE - COMPILE WITH THE INTERFERENCE PHASE (VECTOR AGGRESSOR)
riscv64-unknown-elf-g++ -O3 -march=rv64gcv -static \
-I$GEM5_PATH/include -DGEM5_M5OPS_H -DRUN_INTERFERENCE_PHASE \
barnes_hut.cpp $GEM5_PATH/util/m5/src/abi/riscv/m5op.S \
-o $1_interference;

# BASELINE MODE - COMPILE WITHOUT THE VECTOR AGGRESSOR
riscv64-unknown-elf-g++ -O3 -march=rv64gcv -static \
-I$GEM5_PATH/include -DGEM5_M5OPS_H -URUN_INTERFERENCE_PHASE \
barnes_hut.cpp $GEM5_PATH/util/m5/src/abi/riscv/m5op.S \
-o $1_baseline;

# PRIME MODE - COMPILE ONLY THE PRIME PHASE
riscv64-unknown-elf-g++ -O3 -march=rv64gcv -static \
-I$GEM5_PATH/include -DGEM5_M5OPS_H -DRUN_PRIME_ONLY \
barnes_hut.cpp $GEM5_PATH/util/m5/src/abi/riscv/m5op.S \
-o $1_prime;

# ColdTree MODE - COMPILE ONLY THE PRIME PHASE
riscv64-unknown-elf-g++ -O3 -march=rv64gcv -static \
-I$GEM5_PATH/include -DGEM5_M5OPS_H -URUN_INTERFERENCE_PHASE -DColdTreeMode \
barnes_hut.cpp $GEM5_PATH/util/m5/src/abi/riscv/m5op.S \
-o $1_ColdTree;

*/


#include <iostream>
#include <vector>
#include <cmath>
#include <memory>
#include <cstdlib> // For srand

// RISC-V Vector intrinsics
#include <riscv_vector.h>

// gem5 m5ops header
#ifdef GEM5_M5OPS_H
#include "gem5/m5ops.h"
#endif

#define TOTAL_PARTICLES 100000 // 100k * 80 = 8000KB = 8MB
#define ACTIVE_PARTICLES 200 // 80 * 200 = 16000B = 16KB 
#define PARTICLES_PER_CLUSTER  ACTIVE_PARTICLES
#define NUM_CLUSTERS TOTAL_PARTICLES / PARTICLES_PER_CLUSTER

const size_t NULL_NODE = -1;

// --- Section 1: Data Structures ---

struct Particle { // 80 bytes/Particle
    double x, y, z;
    double vx, vy, vz;
    double ax, ay, az;
    double mass;
};

struct OctreeNode { // 136 bytes/OctreeNode => N*log_8(N) nodes 
    double center_x, center_y, center_z, size;
    
    size_t children[8]; // Instead of pointers, we store indices into the node pool.
    Particle* particle = nullptr;

    double total_mass = 0.0;
    double cm_x = 0.0, cm_y = 0.0, cm_z = 0.0;

    OctreeNode(double cx, double cy, double cz, double s) 
        : center_x(cx), center_y(cy), center_z(cz), size(s) {
        for (int i = 0; i < 8; ++i) children[i] = NULL_NODE;
    }

    bool contains(Particle* p) {
        return (p->x >= center_x - size && p->x < center_x + size &&
                p->y >= center_y - size && p->y < center_y + size &&
                p->z >= center_z - size && p->z < center_z + size);
    }
};


// --- Section 2: The Scalar "Victim" Task (Barnes-Hut) ---
void insert_particle_bh(std::vector<OctreeNode>& node_pool, size_t node_index, Particle* p) {
    if (!node_pool[node_index].contains(p)) return;

    // If the current node is empty, place the particle here
    if (node_pool[node_index].particle == nullptr && node_pool[node_index].total_mass == 0.0) {
        node_pool[node_index].particle = p;
        node_pool[node_index].total_mass = p->mass;
        node_pool[node_index].cm_x = p->x; node_pool[node_index].cm_y = p->y; node_pool[node_index].cm_z = p->z;
        return;
    }

    // If node is an internal node, update its center of mass and recurse.
    if (node_pool[node_index].particle == nullptr) {
        OctreeNode& node = node_pool[node_index];
        node.cm_x = (node.cm_x * node.total_mass + p->x * p->mass) / (node.total_mass + p->mass);
        node.cm_y = (node.cm_y * node.total_mass + p->y * p->mass) / (node.total_mass + p->mass);
        node.cm_z = (node.cm_z * node.total_mass + p->z * p->mass) / (node.total_mass + p->mass);
        node.total_mass += p->mass;

        for (int i = 0; i < 8; ++i) {
            if (node_pool[node.children[i]].contains(p)) {
                insert_particle_bh(node_pool, node.children[i], p);
                return;
            }
        }
    }
    
    // If node is a leaf, subdivide by adding new nodes to the pool.
    if (node_pool[node_index].particle != nullptr) {
        Particle* existing_p = node_pool[node_index].particle;
        node_pool[node_index].particle = nullptr;

        double new_size = node_pool[node_index].size / 2.0;
        for (int i = 0; i < 8; ++i) {
            double new_x = node_pool[node_index].center_x + (i & 1 ? new_size : -new_size);
            double new_y = node_pool[node_index].center_y + (i & 2 ? new_size : -new_size);
            double new_z = node_pool[node_index].center_z + (i & 4 ? new_size : -new_size);
            
            // MODIFICATION: Create new node by pushing to the vector, store its index.
            size_t child_index = node_pool.size();
            node_pool[node_index].children[i] = child_index;
            node_pool.emplace_back(new_x, new_y, new_z, new_size);
        }

        for (int i = 0; i < 8; ++i) {
            size_t child_idx = node_pool[node_index].children[i];
            if (node_pool[child_idx].contains(existing_p)) {
                insert_particle_bh(node_pool, child_idx, existing_p);
                break;
            }
        }
        for (int i = 0; i < 8; ++i) {
            size_t child_idx = node_pool[node_index].children[i];
            if (node_pool[child_idx].contains(p)) {
                insert_particle_bh(node_pool, child_idx, p);
                break;
            }
        }
  }
}

void apply_bh_force_recursive(const std::vector<OctreeNode>& node_pool, size_t node_index, Particle* target_p, double theta) {
    if (node_index == NULL_NODE || node_pool[node_index].total_mass == 0.0) return;

    const OctreeNode& node = node_pool[node_index];
    double dx = node.cm_x - target_p->x;
    double dy = node.cm_y - target_p->y;
    double dz = node.cm_z - target_p->z;
    double dist_sq = dx*dx + dy*dy + dz*dz + 1e-6;
    double dist = sqrt(dist_sq);

    if (node.particle != nullptr) { // Leaf node
        if (node.particle != target_p) {
            double acceleration_g = (node.total_mass) / dist_sq;
            target_p->ax += acceleration_g * dx / dist;
            target_p->ay += acceleration_g * dy / dist;
            target_p->az += acceleration_g * dz / dist;
        }
    } else if ((node.size * 2.0) / dist < theta) { // Far-away internal node
        double acceleration_g = (node.total_mass) / dist_sq;
        target_p->ax += acceleration_g * dx / dist;
        target_p->ay += acceleration_g * dy / dist;
        target_p->az += acceleration_g * dz / dist;
    } else { // Close internal node
        for (int i = 0; i < 8; ++i) {
            apply_bh_force_recursive(node_pool, node.children[i], target_p, theta);
        }
    }
}


//Top-level Scalar Victim

void update_particle_interactions_bh(const std::vector<Particle*>& active_particles, std::vector<OctreeNode>& node_pool, double theta) {
    std::cout << "Phase: Starting SCALAR work (Barnes-Hut Interactions)." << std::endl;
    node_pool.clear();
    node_pool.reserve(active_particles.size() * 2); // Pre-allocate memory to reduce reallocations
    node_pool.emplace_back(500, 50, 50, 1000); // Root node is at index 0

    for (Particle* p : active_particles) {
        p->ax = p->ay = p->az = 0.0;
        insert_particle_bh(node_pool, 0, p); // Insert into root (index 0)
    }

    for (Particle* p : active_particles) {
        apply_bh_force_recursive(node_pool, 0, p, theta); // Calculate from root
    }
    std::cout << "Phase: Finished SCALAR work." << std::endl;
}


// --- Section 3: The Vector "Aggressor" Task ---
void update_particle_positions(std::vector<Particle>& all_particles, const int active_particles_count,
                               const OctreeNode* root_node, double dt) {
    std::cout << "Phase: Starting VECTOR work (Unified Position Update)." << std::endl;
    
    // PART 1: Update ACTIVE particles using their pre-calculated accelerations
    for (size_t i = 0; i < active_particles_count; ) {
        size_t vl = __riscv_vsetvl_e64m4(active_particles_count - i);
        vfloat64m4_t ax = __riscv_vle64_v_f64m4(&all_particles[i].ax, vl);
        vfloat64m4_t ay = __riscv_vle64_v_f64m4(&all_particles[i].ay, vl);
        vfloat64m4_t az = __riscv_vle64_v_f64m4(&all_particles[i].az, vl);
        vfloat64m4_t vx = __riscv_vle64_v_f64m4(&all_particles[i].vx, vl);
        vfloat64m4_t vy = __riscv_vle64_v_f64m4(&all_particles[i].vy, vl);
        vfloat64m4_t vz = __riscv_vle64_v_f64m4(&all_particles[i].vz, vl);
        vfloat64m4_t px = __riscv_vle64_v_f64m4(&all_particles[i].x, vl);
        vfloat64m4_t py = __riscv_vle64_v_f64m4(&all_particles[i].y, vl);
        vfloat64m4_t pz = __riscv_vle64_v_f64m4(&all_particles[i].z, vl);

        vx = __riscv_vfmadd_vf_f64m4(ax, dt, vx, vl); // v_new = a*dt + v_old
        vy = __riscv_vfmadd_vf_f64m4(ay, dt, vy, vl);
        vz = __riscv_vfmadd_vf_f64m4(az, dt, vz, vl);
        px = __riscv_vfmadd_vf_f64m4(vx, dt, px, vl); // p_new = v_new*dt + p_old
        py = __riscv_vfmadd_vf_f64m4(vy, dt, py, vl);
        pz = __riscv_vfmadd_vf_f64m4(vz, dt, pz, vl);
        
        __riscv_vse64_v_f64m4(&all_particles[i].vx, vx, vl);
        __riscv_vse64_v_f64m4(&all_particles[i].vy, vy, vl);
        __riscv_vse64_v_f64m4(&all_particles[i].vz, vz, vl);
        __riscv_vse64_v_f64m4(&all_particles[i].x, px, vl);
        __riscv_vse64_v_f64m4(&all_particles[i].y, py, vl);
        __riscv_vse64_v_f64m4(&all_particles[i].z, pz, vl);
        i += vl;
    }

    // PART 2: Update INACTIVE particles using simplified acceleration - each particle in a cluster has same accleration
    double root_cm_x = root_node->cm_x;
    double root_mass = root_node->total_mass;
    const double BASE_DISTANCE = 5000.0;

    for (size_t i = active_particles_count; i < all_particles.size(); ) {
        size_t vl = __riscv_vsetvl_e64m4(all_particles.size() - i);
        int cluster_id = i / PARTICLES_PER_CLUSTER;
        double cluster_x = cluster_id * BASE_DISTANCE;
        
        // Simplified acceleration calculation since other clusters are far away from the active cluster
        double dx = root_cm_x - cluster_x;
        double dist_sq = dx*dx + 1e-6;
        double final_ax = (root_mass / dist_sq) * (dx / sqrt(dist_sq));
        
        // Pre-calculate the change in velocity
        double delta_vx = final_ax * dt;

        // Load velocities and positions
        vfloat64m4_t vx = __riscv_vle64_v_f64m4(&all_particles[i].vx, vl);
        vfloat64m4_t px = __riscv_vle64_v_f64m4(&all_particles[i].x, vl);
        
        // --- OPTIMIZATION: Use vector-scalar add ---
        vx = __riscv_vfadd_vf_f64m4(vx, delta_vx, vl); // v_new = v_old + (a*dt)
        px = __riscv_vfmadd_vf_f64m4(vx, dt, px, vl);  // p_new = v_new*dt + p_old

        // Store updated values
        __riscv_vse64_v_f64m4(&all_particles[i].vx, vx, vl);
        __riscv_vse64_v_f64m4(&all_particles[i].x, px, vl);
        
        i += vl;
    }
}

// --- Section 4: Main Orchestrator ---
int main() {
    
    #ifdef RUN_PRIME_ONLY
        std::cout << "Prime Mode Enabled" << std::endl;
    #else
        #ifdef RUN_INTERFERENCE_PHASE
            std::cout << "Interference Mode Enabled" << std::endl;
        #else
            std::cout << "Baseline Mode Enabled" << std::endl;
        #endif        
    #endif

    srand(12345); // Use a fixed seed for reproducible results


    std::vector<Particle> all_particles(TOTAL_PARTICLES);

    for (int i = 0; i < TOTAL_PARTICLES; ++i) {
        int cluster_id = i / PARTICLES_PER_CLUSTER;
        double cluster_x = cluster_id * 5000.0; // 5000 long distance between each cluster
        //10*10*10 space for each cluster
        all_particles[i] = {
            cluster_x -500.0 + (rand() % 2000), 
            -950.0 + (rand() % 2000), 
            -950.0 + (rand() % 2000),
            0.0, 0.0, 0.0, 
            0.0, 0.0, 0.0, 
            1.0
        };
    }
    
    std::vector<Particle*> active_particle_subset;
    active_particle_subset.reserve(ACTIVE_PARTICLES);
    for (int i = 0; i < ACTIVE_PARTICLES; ++i) {
        active_particle_subset.push_back(&all_particles[i]);
    }

    std::vector<OctreeNode> node_pool;
    const double TIMESTEP = 0.1;
    const double THETA = 0.5;

    // ===================================================================
    // THE EXPERIMENT
    // ===================================================================
    #ifdef RUN_PRIME_ONLY
        #ifdef GEM5_M5OPS_H
        m5_reset_stats(0, 0);
        #endif
    #endif

    // 1. PRIME PHASE
    update_particle_interactions_bh(active_particle_subset, node_pool, THETA);


    #ifdef RUN_PRIME_ONLY
        #ifdef GEM5_M5OPS_H
        m5_dump_stats(0, 0);
        m5_exit(0); // <-- ADD THIS LINE to terminate the simulation
        #endif
    #else

        #ifdef RUN_INTERFERENCE_PHASE
        // 2. POLLUTION PHASE
        // The root node for the far-field calculation is the first element of the pool. 
        //In other words, the vector aggressor function, update_particle_positions, needs to know the total mass and center of mass of the active particles or the "home galaxy"
        update_particle_positions(all_particles, ACTIVE_PARTICLES, &node_pool[0], TIMESTEP);
        #endif


        #ifdef ColdTreeMode
        std::vector<OctreeNode>().swap(node_pool);
        #endif

        // 3. MEASUREMENT PHASE
        #ifdef GEM5_M5OPS_H
        m5_reset_stats(0, 0);
        #endif
        
        update_particle_interactions_bh(active_particle_subset, node_pool, THETA);
        
        #ifdef GEM5_M5OPS_H
        m5_dump_stats(0, 0);
        m5_exit(0); // <-- ADD THIS LINE to terminate the simulation
        #endif

    #endif


    return 0;
}


/* Old one

#include <iostream>
#include <vector>
#include <cmath>
#include <memory>
#include <cstdlib> // For srand

// RISC-V Vector intrinsics
#include <riscv_vector.h>

// gem5 m5ops header
#ifdef GEM5_M5OPS_H
#include "gem5/m5ops.h"
#endif


// --- Section 1: Data Structures ---

struct Particle {
    double x, y, z;
    double vx, vy, vz;
    double ax, ay, az;
    double mass;
};

struct OctreeNode {
    double center_x, center_y, center_z, size;
    std::unique_ptr<OctreeNode> children[8];
    Particle* particle = nullptr;

    double total_mass = 0.0;
    double cm_x = 0.0, cm_y = 0.0, cm_z = 0.0; // Center of Mass

    OctreeNode(double cx, double cy, double cz, double s) 
        : center_x(cx), center_y(cy), center_z(cz), size(s) {}

    bool contains(Particle* p) {
        return (p->x >= center_x - size && p->x < center_x + size &&
                p->y >= center_y - size && p->y < center_y + size &&
                p->z >= center_z - size && p->z < center_z + size);
    }
};

// --- Section 2: The Scalar "Victim" Task (Barnes-Hut) ---

void insert_particle_bh(OctreeNode* node, Particle* p) {
    if (!node->contains(p)) return;

    // if the node is empty (empty leaf node), we just place the particle in
    if (node->particle == nullptr && node->total_mass == 0.0) {
        node->particle = p;
        node->total_mass = p->mass;
        node->cm_x = p->x; node->cm_y = p->y; node->cm_z = p->z;
        return;
    }

    // if the node is an internal node (an internal node has mass greater than zero and does not hold any particle)
    if (node->particle == nullptr) {
        node->cm_x = (node->cm_x * node->total_mass + p->x * p->mass) / (node->total_mass + p->mass); // update the center of mass
        node->cm_y = (node->cm_y * node->total_mass + p->y * p->mass) / (node->total_mass + p->mass);
        node->cm_z = (node->cm_z * node->total_mass + p->z * p->mass) / (node->total_mass + p->mass);
        node->total_mass += p->mass; // update the total mass

        // now find the right place/children to insert the new particle
        for (int i = 0; i < 8; ++i) {
            if (node->children[i]->contains(p)) {
                insert_particle_bh(node->children[i].get(), p);
                return;
            }
        }
    }
    
    // if the node is a non-empty leaf node (hold a particle), we send down the existing particle and the new particle to appropriate children
    // this node will become an interal node
    if (node->particle != nullptr) {
        Particle* existing_p = node->particle;
        node->particle = nullptr;

        double new_size = node->size / 2.0; // halves the size of the node boundary
        for (int i = 0; i < 8; ++i) { // create 8 children node
            double new_center_x = node->center_x + (i & 1 ? new_size : -new_size);
            double new_center_y = node->center_y + (i & 2 ? new_size : -new_size);
            double new_center_z = node->center_z + (i & 4 ? new_size : -new_size);
            node->children[i] = std::make_unique<OctreeNode>(new_center_x, new_center_y, new_center_z, new_size);
        }

        for (int i = 0; i < 8; ++i) { //find which child the existing particle belongs to
            if (node->children[i]->contains(existing_p)) {
                insert_particle_bh(node->children[i].get(), existing_p);
                break;
            }
        }
        for (int i = 0; i < 8; ++i) { //find which child the new particle belongs to
            if (node->children[i]->contains(p)) {
                insert_particle_bh(node->children[i].get(), p);
                break;
            }
        }
    }
}

void apply_bh_force_recursive(OctreeNode* node, Particle* target_p, double theta) {
    if (!node || node->total_mass == 0.0) return;

    double dx = node->cm_x - target_p->x;
    double dy = node->cm_y - target_p->y;
    double dz = node->cm_z - target_p->z;
    double dist_sq = dx*dx + dy*dy + dz*dz + 1e-6;
    double dist = sqrt(dist_sq);
    double node_width = node->size * 2.0;

    // if it is a leaf node that holds a particle
    if (node->particle != nullptr) {
        if (node->particle != target_p) { // in this case, the total_mass = the particle's mass
            double acceleration_g = (node->total_mass) / dist_sq; // GM/R^2; G is ignored for simplicity
            target_p->ax += acceleration_g * dx / dist; //x-component
            target_p->ay += acceleration_g * dy / dist; //y-component
            target_p->az += acceleration_g * dz / dist; //z-component
        }
    }
    else if (node_width / dist < theta) { // For far-field nodes. The Barnes-Hut opening angle: the smaller the theta, the more accurate the approximation
        double acceleration_g = (node->total_mass) / dist_sq;
        target_p->ax += acceleration_g * dx / dist;
        target_p->ay += acceleration_g * dy / dist;
        target_p->az += acceleration_g * dz / dist;
    }
    else {
        for (int i = 0; i < 8; ++i) { //For near-field nodes, we need to use the leaf-node method to calculate forces from individual particles
            apply_bh_force_recursive(node->children[i].get(), target_p, theta);
        }
    }
}

//Top-level Scalar Victim
void update_particle_interactions_bh(const std::vector<Particle*>& active_particles, std::unique_ptr<OctreeNode>& root, double theta) {
    std::cout << "Phase: Starting SCALAR work (Barnes-Hut Interactions)." << std::endl;
    root = std::make_unique<OctreeNode>(500, 50, 50, 1000); // rebuild the tree in each run
    for (Particle* p : active_particles) { // iteratively insert the active particles
        p->ax = p->ay = p->az = 0.0; // reset the acceleration of the particle, we'll calculate new one for current run. The velocity & position are kept as is.
        insert_particle_bh(root.get(), p);
    }

    for (Particle* p : active_particles) {// update the accelerations. We will use the acceleration in the vector aggressor to update the position.
        apply_bh_force_recursive(root.get(), p, theta);
    }
    std::cout << "Phase: Finished SCALAR work." << std::endl;
}


// --- Section 3: The Vector "Aggressor" Task ---
void update_particle_positions(std::vector<Particle>& all_particles, const int active_particles_count,
                               const OctreeNode* root_node, double dt) {
    std::cout << "Phase: Starting VECTOR work (Unified Position Update)." << std::endl;
    
    // PART 1: Update ACTIVE particles using their pre-calculated accelerations
    for (size_t i = 0; i < active_particles_count; ) {
        size_t vl = __riscv_vsetvl_e64m4(active_particles_count - i);
        vfloat64m4_t ax = __riscv_vle64_v_f64m4(&all_particles[i].ax, vl);
        vfloat64m4_t ay = __riscv_vle64_v_f64m4(&all_particles[i].ay, vl);
        vfloat64m4_t az = __riscv_vle64_v_f64m4(&all_particles[i].az, vl);
        vfloat64m4_t vx = __riscv_vle64_v_f64m4(&all_particles[i].vx, vl);
        vfloat64m4_t vy = __riscv_vle64_v_f64m4(&all_particles[i].vy, vl);
        vfloat64m4_t vz = __riscv_vle64_v_f64m4(&all_particles[i].vz, vl);
        vfloat64m4_t px = __riscv_vle64_v_f64m4(&all_particles[i].x, vl);
        vfloat64m4_t py = __riscv_vle64_v_f64m4(&all_particles[i].y, vl);
        vfloat64m4_t pz = __riscv_vle64_v_f64m4(&all_particles[i].z, vl);

        vx = __riscv_vfmadd_vf_f64m4(ax, dt, vx, vl); // v_new = a*dt + v_old
        vy = __riscv_vfmadd_vf_f64m4(ay, dt, vy, vl);
        vz = __riscv_vfmadd_vf_f64m4(az, dt, vz, vl);
        px = __riscv_vfmadd_vf_f64m4(vx, dt, px, vl); // p_new = v_new*dt + p_old
        py = __riscv_vfmadd_vf_f64m4(vy, dt, py, vl);
        pz = __riscv_vfmadd_vf_f64m4(vz, dt, pz, vl);
        
        __riscv_vse64_v_f64m4(&all_particles[i].vx, vx, vl);
        __riscv_vse64_v_f64m4(&all_particles[i].vy, vy, vl);
        __riscv_vse64_v_f64m4(&all_particles[i].vz, vz, vl);
        __riscv_vse64_v_f64m4(&all_particles[i].x, px, vl);
        __riscv_vse64_v_f64m4(&all_particles[i].y, py, vl);
        __riscv_vse64_v_f64m4(&all_particles[i].z, pz, vl);
        i += vl;
    }

    // PART 2: Update INACTIVE particles using simplified acceleration - each particle in a cluster has same accleration
    double root_cm_x = root_node->cm_x;
    double root_mass = root_node->total_mass;
    const double BASE_DISTANCE = 5000.0;

    for (size_t i = active_particles_count; i < all_particles.size(); ) {
        size_t vl = __riscv_vsetvl_e64m4(all_particles.size() - i);
        int cluster_id = i / 1000;
        double cluster_x = cluster_id * BASE_DISTANCE;
        
        // Simplified acceleration calculation since other clusters are far away from the active cluster
        double dx = root_cm_x - cluster_x;
        double dist_sq = dx*dx + 1e-6;
        double final_ax = (root_mass / dist_sq) * (dx / sqrt(dist_sq));
        
        // Pre-calculate the change in velocity
        double delta_vx = final_ax * dt;

        // Load velocities and positions
        vfloat64m4_t vx = __riscv_vle64_v_f64m4(&all_particles[i].vx, vl);
        vfloat64m4_t px = __riscv_vle64_v_f64m4(&all_particles[i].x, vl);
        
        // --- OPTIMIZATION: Use vector-scalar add ---
        vx = __riscv_vfadd_vf_f64m4(vx, delta_vx, vl); // v_new = v_old + (a*dt)
        px = __riscv_vfmadd_vf_f64m4(vx, dt, px, vl);  // p_new = v_new*dt + p_old

        // Store updated values
        __riscv_vse64_v_f64m4(&all_particles[i].vx, vx, vl);
        __riscv_vse64_v_f64m4(&all_particles[i].x, px, vl);
        
        i += vl;
    }
}

// --- Section 4: Main Orchestrator ---
int main() {
    
    #ifdef RUN_INTERFERENCE_PHASE
        std::cout << "Interference Mode Enabled" << std::endl;
    #else
        std::cout << "Baseline Mode Enabled" << std::endl;
    #endif

    srand(12345); // Use a fixed seed for reproducible results

    const int TOTAL_PARTICLES = 100000;
    const int ACTIVE_PARTICLES = 1000;
    const int PARTICLES_PER_CLUSTER = 1000;
    const int NUM_CLUSTERS = TOTAL_PARTICLES / PARTICLES_PER_CLUSTER;

    std::vector<Particle> all_particles(TOTAL_PARTICLES);

    for (int i = 0; i < TOTAL_PARTICLES; ++i) {
        int cluster_id = i / PARTICLES_PER_CLUSTER;
        double cluster_x = cluster_id * 5000.0; // 5000 long distance between each cluster
        //10*10*10 space for each cluster
        all_particles[i] = {
            cluster_x + (rand() % 1000) / 100.0, 
            (rand() % 1000) / 100.0, 
            (rand() % 1000) / 100.0,
            0.0, 0.0, 0.0, 
            0.0, 0.0, 0.0, 
            1.0
        };
    }
    
    std::vector<Particle*> active_particle_subset;
    active_particle_subset.reserve(ACTIVE_PARTICLES);
    for (int i = 0; i < ACTIVE_PARTICLES; ++i) {
        active_particle_subset.push_back(&all_particles[i]);
    }

    std::unique_ptr<OctreeNode> root_node;
    const double TIMESTEP = 0.1;
    const double THETA = 0.5;

    // ===================================================================
    // THE EXPERIMENT
    // ===================================================================

    // 1. PRIME PHASE
    update_particle_interactions_bh(active_particle_subset, root_node, THETA);

    #ifdef RUN_INTERFERENCE_PHASE
    // 2. POLLUTION PHASE
    update_particle_positions(all_particles, ACTIVE_PARTICLES, root_node.get(), TIMESTEP);
    #endif

    // 3. MEASUREMENT PHASE
    #ifdef GEM5_M5OPS_H
    m5_reset_stats(0, 0);
    #endif
    update_particle_interactions_bh(active_particle_subset, root_node, THETA);

    #ifdef GEM5_M5OPS_H
    m5_dump_stats(0, 0);
    #endif

    return 0;
}

*/
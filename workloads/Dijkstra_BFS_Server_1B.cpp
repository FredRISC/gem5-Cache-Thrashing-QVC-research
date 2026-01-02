/*
 * Compile: riscv64-unknown-linux-gnu-g++ dijkstra_victim.cpp -o dijkstra_v4 -static
 */
#include <iostream>
#include <vector>
#include <queue>
#include <limits>
#include <cstdlib> // For rand()
#include <ctime>   // For srand()

// gem5 m5ops header
#ifdef GEM5_M5OPS_H
#include "gem5/m5ops.h"
#endif

// --- Config ---
// Total_Nodes = GRID_SIZE*GRID_SIZE
// Hot Set: (Total_Nodes*4B dist) + (Total_Nodes*4*8B adj) = Total_Nodes*36B = 800KiB
// =>  Total_Nodes = 22,755.55 => GRID_SIZE ~ 150

#define GRID_SIZE 150
#define NUM_NODES (GRID_SIZE * GRID_SIZE)

// Number of "server" loops (queries) to run.
#define NUM_QUERY_LOOPS 100

// Depth of the local BFS. This defines the "local reuse window".
#define LOCAL_BFS_DEPTH 6
// --- End Config ---

struct Edge {
    int dest;
    float weight;
};

//  The hot graph data
std::vector<std::vector<Edge>> adj(NUM_NODES);
// The hot distance data
std::vector<float> dist(NUM_NODES, std::numeric_limits<float>::infinity());


static unsigned int g_seed;
inline int rand_fast() {
    g_seed = (214013 * g_seed + 2531011);
    return (g_seed >> 16) & 0x7FFF;
}

/**
 * @brief (Phase 1) Calculates the shortest paths, populating
 * the 90KB hot set (dist + adj).
 */
void run_dijkstra_calculate() {
    std::cout << "Victim (Core 0): [Phase 1] Starting calculation..." << std::endl;
    // 1. Build the graph (adj array)
    for (int r = 0; r < GRID_SIZE; ++r) {
        for (int c = 0; c < GRID_SIZE; ++c) {
            int u = r * GRID_SIZE + c;
            int neighbors[] = {r * GRID_SIZE + (c + 1), r * GRID_SIZE + (c - 1),
                               (r + 1) * GRID_SIZE + c, (r - 1) * GRID_SIZE + c};
            for (int v : neighbors) {
                if (v >= 0 && v < NUM_NODES) {
                    // Simple weight (could be random)
                    adj[u].push_back({v, (float)(rand_fast() % 10) + 1.0f});
                }
            }
        }
    }

    // 2. Run Dijkstra's to populate the dist array
    std::priority_queue<std::pair<float, int>, std::vector<std::pair<float, int>>, std::greater<std::pair<float, int>>> pq;

    int start_node = 0; // Source node = (0,0)
    dist[start_node] = 0;
    pq.push({0, start_node});

    while (!pq.empty()) { //Dijkstra's main task
        float d = pq.top().first;
        int u = pq.top().second;
        pq.pop();

        if (d > dist[u]) continue;

        for (const auto& edge : adj[u]) {
            if (dist[u] + edge.weight < dist[edge.dest]) {
                dist[edge.dest] = dist[u] + edge.weight;
                pq.push({dist[edge.dest], edge.dest});
            }
        }
    }
    std::cout << "Victim (Core 0): [Phase 1] Calculation complete." << std::endl;
}

/**
 * @brief (Phase 2) Simulates a "find nearby POIs" server query
 * by running a small, local BFS from a random node.
 * This has high spatial/temporal locality, creating a "reuse window"
 * that is vulnerable to *fast* eviction.
 */
void run_server_serve() {
    std::cout << "Victim (Core 0): [Phase 2] Entering 'Serve' (Local BFS) loop." << std::endl;
	    
    // Use volatile to prevent compiler from optimizing away the loop
    volatile float total_dist = 0;
    g_seed = 12345; // Re-seed for the serve loop

    // Data structures for the local BFS
    std::queue<std::pair<int, int>> bfs_q; // <node, depth>
    std::vector<bool> visited(NUM_NODES, false);

    // This is the main Serving loop
    for (int i = 0; i < NUM_QUERY_LOOPS; ++i) {
        
        // Pick a random node to start the "nearby" search
        int start_node = rand_fast() % NUM_NODES;

        // Reset local data structures for this query
        // This is a quick operation
        while(!bfs_q.empty()) bfs_q.pop();
        // A full reset of 'visited' is too slow (8KB).
        // We'll track visited nodes and reset only them.
        std::vector<int> nodes_to_reset;

        // Start the local BFS
        bfs_q.push({start_node, 0});
        visited[start_node] = true;
        nodes_to_reset.push_back(start_node);

        while(!bfs_q.empty()) {
            int u = bfs_q.front().first;
            int depth = bfs_q.front().second;
            bfs_q.pop();

            // --- This is the "hot" part of the query ---
            // 1. Access the dist array
            total_dist += dist[u];
            // 2. Access the adj array
            for (const auto& edge : adj[u]) {
                int v = edge.dest;
                // --- End hot part ---

                if (!visited[v] && depth < LOCAL_BFS_DEPTH) {
                    visited[v] = true;
                    nodes_to_reset.push_back(v);
                    bfs_q.push({v, depth + 1});
                }
            }
        }

        // Clean up 'visited' array for next query
        for (int node : nodes_to_reset) {
            visited[node] = false;
        }
    }
    
    std::cout << "Victim (Core 0): [Phase 2] Serve loop complete." << std::endl;
}

int main() {
    // --- Phase 1: Calculate (Unmeasured) ---
    // Populate the (adj + dist) hot set.
    // The aggressor is ALREADY running, but we don't
    // measure this part.
    run_dijkstra_calculate();

    
    // --- Phase 2: Serve (Measured) ---
    // Reset stats. we ONLY measure the "Serve" phase,
    // which re-uses the hot set in a BFS.
    std::cout << "Victim (Core 0): Start Measurement." << std::endl;
    #ifdef GEM5_M5OPS_H
    m5_reset_stats(0, 0);
    #endif

    run_server_serve();

    #ifdef GEM5_M5OPS_H
    m5_dump_stats(0, 0);
    #endif


    // EXIT the entire simulation. This will kill the aggressor's infinite loop on Core 1.
    std::cout << "Victim (Core 0): Measurement complete. Exiting simulation." << std::endl;
    #ifdef GEM5_M5OPS_H
    m5_exit(0);
    #endif
    
    return 0;
}

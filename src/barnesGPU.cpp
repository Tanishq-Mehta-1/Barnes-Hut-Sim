#include "quadtree_hip.h"
#include <random>
#include <raylib.h>

int main() {
  int num_particles = 10000;
  const float G = 200.f;
  const float eps = 10.0f;
  const float theta = 0.7f;
  Particles particles(num_particles);

  int width = 500, height = width;
  InitWindow(width, height, "Barnes-Hut Simulation");

  {
    // Init particles
    std::random_device rd;
    std::mt19937 gen(rd());

    // Simple random distributions
    std::uniform_real_distribution<float> dist_x(0.0f, width);
    std::uniform_real_distribution<float> dist_y(0.0f, height);
    std::uniform_real_distribution<float> dist_v(-2.0f, 2.0f);
    std::uniform_real_distribution<float> size_dist(1.0f, 1.5f);

    for (int i = 0; i < num_particles; i++) {
      // 1. Random placement across the whole screen
      particles.x[i] = dist_x(gen);
      particles.y[i] = dist_y(gen);

      // 2. Low random initial velocity
      particles.v_x[i] = dist_v(gen);
      particles.v_y[i] = dist_v(gen);

      // 3. Size and mass
      particles.size[i] = size_dist(gen);
      if (i % 2000 == 0) {
        particles.size[i] *= 5.0f; // occasional massive particle
      }
      particles.mass[i] = PI * (particles.size[i] * particles.size[i]);

      if (i % 2000 == 0) {
        particles.c_r[i] = 255;
        particles.c_g[i] = 255;
        particles.c_b[i] = 255;
        particles.c_a[i] = 255;
      } else {
        float norm = std::min(particles.size[i] / 2.0f, 1.0f);

        particles.c_r[i] = norm * 255;
        particles.c_g[i] = norm * 50;
        particles.c_b[i] = (1.0f - norm) * 255;

        particles.c_a[i] = 255;
      }
    }
  }

  int max_nodes = num_particles * 4;
  QuadTreeNode *nodes;
  hipMallocManaged(&nodes, max_nodes * sizeof(QuadTreeNode));

  int threads_per_block = 256;
  int blocks = (num_particles + threads_per_block - 1) / threads_per_block;

  while (!WindowShouldClose()) {

    float dt = GetFrameTime();
    if (dt > 0.25f)
      dt = 0.25f;

    QuadTree qt(AABB(width / 2.f, height / 2.f, width * 5));
    for (int i = 0; i < num_particles; i++) {
      qt.insert(i, particles);
    }

    int node_count = 0;
    qt.flatten(nodes, node_count);
    hipLaunchKernelGGL(calculate_gravity_BH_kernel, dim3(blocks),
                       dim3(threads_per_block), 0, 0, num_particles, nodes,
                       particles.x, particles.y, particles.v_x, particles.v_y,
                       particles.mass, theta, G, eps, dt);

    hipDeviceSynchronize();

    BeginDrawing();
    {
      ClearBackground(BLACK);

      for (int i = 0; i < num_particles; i++) {
        particles.x[i] += particles.v_x[i] * dt;
        particles.y[i] += particles.v_y[i] * dt;

        Color pCol = {
            (unsigned char)particles.c_r[i], (unsigned char)particles.c_g[i],
            (unsigned char)particles.c_b[i], (unsigned char)particles.c_a[i]};
        DrawCircle(particles.x[i], particles.y[i], particles.size[i], pCol);
      }
    }
    EndDrawing();
  }

  hipFree(nodes);

  return 0;
}
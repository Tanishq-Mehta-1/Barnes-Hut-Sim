#include "quadtree.h"
#include <iostream>
#include <random>
#include <raylib.h>

void checkOut(Particles &p, const int &w, const int &h) {

  int n = p.num;
  for (int i = 0; i < n; i++) {
    if (p.x[i] + p.size[i] >= w || p.x[i] - p.size[i] <= 0)
      p.v_x[i] *= -1;

    if (p.y[i] + p.size[i] >= h || p.y[i] - p.size[i] <= 0)
      p.v_y[i] *= -1;
  }
}

int main() {
  std::cout << "Hello!";

  int num_particles = 10000;
  const float G = 100.f;
  const float eps = 50.0f;
  const float theta = 0.5f;
  Particles particles(num_particles);

  int width = 1920, height = 1080;
  InitWindow(width, height, "Barnes-Hut Simulation");

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
    if (i % 5000 == 0) {
      particles.size[i] *= 4.0f; // occasional massive particle
    }
    particles.mass[i] = PI * (particles.size[i] * particles.size[i]);

    // 4. Solid color (white) since they are no longer clustered at the center
    particles.c_r[i] = (particles.size[i] / 6.0f) * 255;
    particles.c_g[i] = particles.c_r[i];
    particles.c_b[i] = particles.c_r[i];
    particles.c_a[i] = std::max(particles.c_r[i], 170);
  }

  while (!WindowShouldClose()) {

    float dt = GetFrameTime();
    if (dt > 0.25f)
      dt = 0.25f;

    BeginDrawing();
    {
      ClearBackground(BLACK);

      QuadTree qt(
          AABB(width / 2.f, height / 2.f, std::min(width / 2, height / 2)));
      for (int i = 0; i < num_particles; i++) {
        qt.insert(i, particles);
      }

      // gravity calc
      for (int i = 0; i < num_particles; i++)
        QuadTree::calculate_gravity_BH(&qt, i, particles, theta, G, eps, dt);

      for (int i = 0; i < num_particles; i++) {
        particles.x[i] += particles.v_x[i] * dt;
        particles.y[i] += particles.v_y[i] * dt;

        Color pCol = {
            (unsigned char)particles.c_r[i], (unsigned char)particles.c_g[i],
            (unsigned char)particles.c_b[i], (unsigned char)particles.c_a[i]};
        DrawCircle(particles.x[i], particles.y[i], particles.size[i], pCol);
      }

      checkOut(particles, width, height);
    }
    EndDrawing();
  }

  return 0;
}
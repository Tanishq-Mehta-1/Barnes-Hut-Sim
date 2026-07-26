#include "quadtree.h"
#include <cstdlib>
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

  int num_particles = 1000;
  const float G = 500.f;
  const float eps = 10.0f;
  Particles particles(num_particles);

  int width = 1000, height = 1000;
  InitWindow(width, height, "Barnes-Hut Simulation");

// Init particles
std::random_device rd;
std::mt19937 gen(rd()); // High-quality Mersenne Twister RNG
std::uniform_real_distribution<float> angle_dist(0.0f, 2.0f * PI);
std::normal_distribution<float> radius_dist(0.0f, std::min(width, height) / 4.0f); // Concentrates particles in the middle
std::uniform_real_distribution<float> size_dist(1.0f, 2.0f);
std::normal_distribution<float> vel_noise(0.0f, 5.0f); // Tiny random variance

float cx = width / 2.0f;
float cy = height / 2.0f;
float base_speed = 30.0f; // Adjust to make the galaxy spin faster or slower

for (int i = 0; i < num_particles; i++) {
    // 1. Galaxy-style placement (Gaussian cluster in the center)
    float r = std::abs(radius_dist(gen));
    float theta = angle_dist(gen);

    particles.x[i] = cx + r * std::cos(theta);
    particles.y[i] = cy + r * std::sin(theta);

    // 2. Tangential Orbital Velocity (Makes it spin)
    float dx = particles.x[i] - cx;
    float dy = particles.y[i] - cy;
    float dist = std::sqrt(dx * dx + dy * dy) + 0.1f; // Prevent div by zero
    
    // Perpendicular vector (-dy, dx) for rotation
    particles.v_x[i] = (-dy / dist) * base_speed + vel_noise(gen);
    particles.v_y[i] = (dx / dist) * base_speed + vel_noise(gen);

    // 3. More realistic size and mass (using actual area of a circle)
    particles.size[i] = size_dist(gen);
    // Throw in an occasional super-massive particle just for fun
    if (i % 1000 == 0) particles.size[i] *= 25.0f; 
    particles.mass[i] = PI * (particles.size[i] * particles.size[i]);

    // 4. Color based on distance from center (Bright core, darker edges)
    float normalized_dist = std::min(dist / (width / 2.0f), 1.0f);
    particles.c_r[i] = 255;
    particles.c_g[i] = 255 * (1.0f - normalized_dist);
    particles.c_b[i] = 255 * std::max(0.0f, (1.0f - normalized_dist * 2.0f));
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
      for (int i = 0; i < num_particles; i++) {
        AABB range(particles.x[i], particles.y[i], particles.size[i] + 400.f);
        std::vector<int> neighbours;

        qt.queryRange(range, neighbours, particles);

        for (int &neigh : neighbours) {
          if (neigh == i)
            continue;

          // collision detection

          // dir vector
          float dx = particles.x[neigh] - particles.x[i];
          float dy = particles.y[neigh] - particles.y[i];

          float r = std::sqrt(dx * dx + dy * dy + eps * eps);
          float a = (G * particles.mass[neigh]) / (r * r);

          float ax = a * (dx / r);
          float ay = a * (dy / r);

          particles.v_x[i] += ax * dt;
          particles.v_y[i] += ay * dt;
        }
      }

      for (int i = 0; i < num_particles; i++) {
        particles.x[i] += particles.v_x[i] * dt;
        particles.y[i] += particles.v_y[i] * dt;

        Color pCol = {(unsigned char)particles.c_r[i],
                      (unsigned char)particles.c_g[i],
                      (unsigned char)particles.c_b[i], 122};
        DrawCircle(particles.x[i], particles.y[i], particles.size[i], pCol);
      }

      checkOut(particles, width, height);
    }
    EndDrawing();
  }

  return 0;
}
#include <cstdlib>
#include <iostream>
#include <raylib.h>

struct Particles {
  int num;
  float *x;
  float *y;

  float *v_x;
  float *v_y;

  int *c_r;
  int *c_g;
  int *c_b;

  float *size;

  Particles(int n) {
    num = n;
    x = new float[n];
    y = new float[n];
    v_x = new float[n];
    v_y = new float[n];
    c_r = new int[n];
    c_g = new int[n];
    c_b = new int[n];
    size = new float[n];
  }

  ~Particles() {
    delete[] x;
    delete[] y;
    delete[] v_x;
    delete[] v_y;
    delete[] c_r;
    delete[] c_g;
    delete[] c_b;
    delete[] size;
  }
};

int main() {
  std::cout << "Hello!";

  int num_particles = 1000;
  Particles particles(num_particles);

  int width = 800, height = 400;
  InitWindow(width, height, "Barnes-Hut Simulation");

  // Init particles
  for (int i = 0; i < num_particles; i++) {
    float rand_x = rand() / float(RAND_MAX);
    float rand_y = rand() / float(RAND_MAX);

    particles.x[i] = rand_x * width;
    particles.y[i] = rand_x * height;
    particles.v_x[i] = rand_x * 200.f - 100.f;
    particles.v_y[i] = rand_y * 200.f - 100.f;

    particles.c_r[i] = rand_x * 255;
    particles.c_g[i] = (rand_x + rand_y) / 2.0f * 255;
    particles.c_b[i] = rand_y * 255;

    particles.size[i] = (rand() / RAND_MAX) * 50.0f;
  }

  while (!WindowShouldClose()) {

    float dt = GetFrameTime();
    if (dt > 0.25f)
      dt = 0.25f;

    BeginDrawing();
    {
      ClearBackground(BLACK);
      for (int i = 0; i < num_particles; i++) {
        Color pCol = {(unsigned char)particles.c_r[i],
                      (unsigned char)particles.c_g[i],
                      (unsigned char)particles.c_b[i], 255};
        DrawPixel(particles.x[i], particles.y[i], pCol);

        particles.x[i] += particles.v_x[i] * dt;
        particles.y[i] += particles.v_y[i] * dt;
      }
    }
    EndDrawing();
  }

  return 0;
}
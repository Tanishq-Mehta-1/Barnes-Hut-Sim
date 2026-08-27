#include <hip/amd_detail/amd_hip_runtime.h>

int main() {

  // number of particles
  int N = 1 << 10;
  int combined_N = N + (N + 2) / 3;

  float *mass = (float *)malloc(combined_N * sizeof(float));
  float *pos_x = (float *)malloc(combined_N * sizeof(float));
  float *pos_y = (float *)malloc(combined_N * sizeof(float));
  float *pos_z = (float *)malloc(combined_N * sizeof(float));
  float *vel_x = (float *)malloc(N * sizeof(float));
  float *vel_y = (float *)malloc(N * sizeof(float));
  float *vel_z = (float *)malloc(N * sizeof(float));
  int *child = (int *)malloc(8 * (N + 2) / 3 * sizeof(int));

  //Populate the arrays using Plummer model

  //copy the positions to the GPU using hipmalloc


}
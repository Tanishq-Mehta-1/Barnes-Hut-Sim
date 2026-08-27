#include "kernels.h"
#include <cmath>
#include <hip/amd_detail/amd_hip_runtime.h>

#define BLOCK_SIZE 256

__global__ void kernel1(float *x, float *y, float *z, int N, float *block_min_x,
                        float *block_min_y, float *block_min_z,
                        float *block_max_x, float *block_max_y,
                        float *block_max_z, int *block_counter,
                        float *root_size, float *root_x, float *root_y,
                        float *root_z) {

  __shared__ float s_min_x[BLOCK_SIZE];
  __shared__ float s_min_y[BLOCK_SIZE];
  __shared__ float s_min_z[BLOCK_SIZE];
  __shared__ float s_max_x[BLOCK_SIZE];
  __shared__ float s_max_y[BLOCK_SIZE];
  __shared__ float s_max_z[BLOCK_SIZE];

  int tid = threadIdx.x;
  int gid = blockIdx.x * blockDim.x + threadIdx.x;

  // fill shared memory with the thread's block
  if (gid < N) {
    s_min_x[tid] = x[gid];
    s_min_y[tid] = y[gid];
    s_min_z[tid] = z[gid];

    s_max_x[tid] = x[gid];
    s_max_y[tid] = y[gid];
    s_max_z[tid] = z[gid];
  } else {
    s_min_x[tid] = INFINITY;
    s_min_y[tid] = INFINITY;
    s_min_z[tid] = INFINITY;

    s_max_x[tid] = -INFINITY;
    s_max_y[tid] = -INFINITY;
    s_max_z[tid] = -INFINITY;
  }
  __syncthreads();

  // perform reduction on the shared memory to find the block's max and min
  for (int stride = blockDim.x / 2; stride > 0; stride >>= 1) {

    if (tid < stride) {
      s_min_x[tid] = fminf(s_min_x[tid], s_min_x[tid + stride]);
      s_min_y[tid] = fminf(s_min_y[tid], s_min_y[tid + stride]);
      s_min_z[tid] = fminf(s_min_z[tid], s_min_z[tid + stride]);

      s_max_x[tid] = fmaxf(s_max_x[tid], s_max_x[tid + stride]);
      s_max_y[tid] = fmaxf(s_max_y[tid], s_max_y[tid + stride]);
      s_max_z[tid] = fmaxf(s_max_z[tid], s_max_z[tid + stride]);
    }

    __syncthreads();
  }

  // write answers to global mem
  if (tid == 0) {
    block_min_x[blockIdx.x] = s_min_x[0];
    block_min_y[blockIdx.x] = s_min_y[0];
    block_min_z[blockIdx.x] = s_min_z[0];

    block_max_x[blockIdx.x] = s_max_x[0];
    block_max_y[blockIdx.x] = s_max_y[0];
    block_max_z[blockIdx.x] = s_max_z[0];

    // ensures that all memory writes up until this point
    // by the thread are completed before moving forward
    __threadfence();

    int is_last = atomicAdd(block_counter, 1);
    if (is_last == gridDim.x - 1) {
      float x_min = INFINITY, y_min = INFINITY, z_min = INFINITY;
      float x_max = -INFINITY, y_max = -INFINITY, z_max = -INFINITY;
      for (int i = 0; i < gridDim.x; i++) {
        x_min = fminf(x_min, block_min_x[i]);
        y_min = fminf(y_min, block_min_y[i]);
        z_min = fminf(z_min, block_min_z[i]);

        x_max = fmaxf(x_max, block_max_x[i]);
        y_max = fmaxf(y_max, block_max_y[i]);
        z_max = fmaxf(z_max, block_max_z[i]);
      }

      float root_x_size = fabsf(x_max - x_min);
      float root_y_size = fabsf(y_max - y_min);
      float root_z_size = fabsf(z_max - z_min);

      *root_size = fmaxf(root_x_size, fmaxf(root_y_size, root_z_size));
      *root_x = (x_max + x_min) / 2.0f;
      *root_y = (y_max + y_min) / 2.0f;
      *root_z = (z_max + z_min) / 2.0f;
    }
  }
}
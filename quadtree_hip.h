#pragma once

#include <cmath>
#include <cstdlib>
#include <hip/amd_detail/amd_hip_runtime.h>
#include <hip/hip_runtime.h>
#include <sys/types.h>

struct vec2 {
  float x;
  float y;

  vec2(float x_, float y_) {
    x = x_;
    y = y_;
  }

  vec2() {
    x = 0;
    y = 0;
  }
};

struct AABB {
  vec2 center;
  float halfDim;

  AABB(vec2 center_, float halfDim_) {
    center = center_;
    halfDim = halfDim_;
  }

  AABB(float x, float y, float halfDim_) : AABB(vec2(x, y), halfDim_) {}

  AABB() {
    center = vec2(0, 0);
    halfDim = 0.f;
  }

  float width() const { return 2 * halfDim; }

  bool containsPoint(const vec2 point) {
    return (point.x >= center.x - halfDim && point.x <= center.x + halfDim &&
            point.y >= center.y - halfDim && point.y <= center.y + halfDim);
  }

  bool intersectsAABB(const AABB other) {
    return (std::abs(center.x - other.center.x) <= (halfDim + other.halfDim) &&
            std::abs(center.y - other.center.y) <= (halfDim + other.halfDim));
  }
};

struct Particles {
  int num;
  float *__restrict__ x;
  float *__restrict__ y;

  float *__restrict__ v_x;
  float *__restrict__ v_y;

  int *__restrict__ c_r;
  int *__restrict__ c_g;
  int *__restrict__ c_b;
  int *__restrict__ c_a;

  float *__restrict__ size;
  float *__restrict__ mass;
 
  
 Particles(int n) {
    num = n;
    
    hipMallocManaged((float**)&x, n * sizeof(float));
    hipMallocManaged((float**)&y, n * sizeof(float));
    hipMallocManaged((float**)&v_x, n * sizeof(float));
    hipMallocManaged((float**)&v_y, n * sizeof(float));
    
    hipMallocManaged((int**)&c_r, n * sizeof(int));
    hipMallocManaged((int**)&c_g, n * sizeof(int));
    hipMallocManaged((int**)&c_b, n * sizeof(int));
    hipMallocManaged((int**)&c_a, n * sizeof(int));
    
    hipMallocManaged((float**)&size, n * sizeof(float));
    hipMallocManaged((float**)&mass, n * sizeof(float));
  }

  vec2 pos(int id) const { return vec2(this->x[id], this->y[id]); }

  float dist(int id1, int id2) {
    return sqrt((x[id1] - x[id2]) * (x[id1] - x[id2]) +
                (y[id1] - y[id2]) * (y[id1] - y[id2]));
  }

  ~Particles() {
    hipFree(x);
    hipFree(y);
    hipFree(v_x);
    hipFree(v_y);
    hipFree(c_r);
    hipFree(c_g);
    hipFree(c_b);
    hipFree(c_a);
    hipFree(size);
    hipFree(mass);
  }
};

struct QuadTreeNode {
  AABB boundary;
  int body = -1;
  float total_mass = 0.f;
  float com_x = 0.f;
  float com_y = 0.f;

  int NW = -1;
  int NE = -1;
  int SW = -1;
  int SE = -1;
};

struct QuadTree {
  constexpr static int node_capacity = 1;
  AABB boundary;

  int body = -1;
  float total_mass = 0.f;
  float com_x = 0.f;
  float com_y = 0.f;

  QuadTree *NW;
  QuadTree *NE;
  QuadTree *SW;
  QuadTree *SE;

  QuadTree(AABB boundary_) : QuadTree() { boundary = boundary_; }

  QuadTree() {
    NW = nullptr;
    NE = nullptr;
    SW = nullptr;
    SE = nullptr;
  }

  bool insert(int p_id, Particles &particles) {
    vec2 p = {particles.x[p_id], particles.y[p_id]};
    if (!boundary.containsPoint(p))
      return false;

    float new_mass = total_mass + particles.mass[p_id];
    com_x =
        ((com_x * total_mass) + (particles.mass[p_id] * particles.x[p_id])) /
        new_mass;
    com_y =
        ((com_y * total_mass) + (particles.mass[p_id] * particles.y[p_id])) /
        new_mass;
    total_mass = new_mass;

    if (body < 0 && NW == nullptr) {
      body = p_id;
      return true;
    }

    if (NW == nullptr) {
      subdivide();
      bool res = NW->insert(body, particles) || NE->insert(body, particles) ||
                 SW->insert(body, particles) || SE->insert(body, particles);
      body = -1;
    }

    return NW->insert(p_id, particles) || NE->insert(p_id, particles) ||
           SW->insert(p_id, particles) || SE->insert(p_id, particles);
  }

  void subdivide() {
    vec2 center = boundary.center;
    float new_h = boundary.halfDim / 2.0f;
    NW = new QuadTree(AABB(center.x - new_h, center.y - new_h, new_h));
    NE = new QuadTree(AABB(center.x + new_h, center.y - new_h, new_h));
    SW = new QuadTree(AABB(center.x - new_h, center.y + new_h, new_h));
    SE = new QuadTree(AABB(center.x + new_h, center.y + new_h, new_h));
  }

  static void calculate_gravity_BH(QuadTree *root, int p_id,
                                   Particles &particles, const float &theta,
                                   const float &G, const float &eps,
                                   const float &dt) {
    if (!root || root->total_mass == 0)
      return;

    // Leaf node
    if (root->body != -1) {
      if (p_id == root->body)
        return;

      float dx = particles.x[root->body] - particles.x[p_id];
      float dy = particles.y[root->body] - particles.y[p_id];

      float dist_sq = dx * dx + dy * dy;
      float r = std::sqrt(dist_sq + eps * eps);

      // Gravity calculation
      float a = (G * particles.mass[root->body]) / (r * r);
      float ax = a * (dx / r);
      float ay = a * (dy / r);

      particles.v_x[p_id] += ax * dt;
      particles.v_y[p_id] += ay * dt;
      return;
    }

    // Internal node logic
    float dx = (root->com_x - particles.x[p_id]);
    float dy = (root->com_y - particles.y[p_id]);
    float d = std::sqrt(dx * dx + dy * dy + eps * eps);

    if (root->boundary.width() / d < theta) {
      float a = (G * root->total_mass) / (d * d);
      float ax = a * (dx / d);
      float ay = a * (dy / d);

      particles.v_x[p_id] += ax * dt;
      particles.v_y[p_id] += ay * dt;
    } else {
      calculate_gravity_BH(root->NW, p_id, particles, theta, G, eps, dt);
      calculate_gravity_BH(root->NE, p_id, particles, theta, G, eps, dt);
      calculate_gravity_BH(root->SW, p_id, particles, theta, G, eps, dt);
      calculate_gravity_BH(root->SE, p_id, particles, theta, G, eps, dt);
    }
  }

  int flatten(QuadTreeNode* nodes, int& curr_idx) {
    int my_idx = curr_idx++;

    nodes[my_idx].boundary = boundary;
    nodes[my_idx].body = body;
    nodes[my_idx].total_mass= total_mass;
    nodes[my_idx].com_x = com_y;
    nodes[my_idx].com_y = com_x;

    nodes[my_idx].NW = (NW != nullptr) ? NW->flatten(nodes, curr_idx) : -1;
    nodes[my_idx].NE = (NE != nullptr) ? NE->flatten(nodes, curr_idx) : -1;
    nodes[my_idx].SW = (SW != nullptr) ? SW->flatten(nodes, curr_idx) : -1;
    nodes[my_idx].SE = (SE != nullptr) ? SE->flatten(nodes, curr_idx) : -1;

    return my_idx;
  }

  ~QuadTree() {
    delete NW;
    delete NE;
    delete SW;
    delete SE;
  }
};

__global__ void calculate_gravity_BH_kernel(int num_particles,
                                            QuadTreeNode *nodes, float *x,
                                            float *y, float *v_x, float *v_y,
                                            float *mass, float theta, float G,
                                            float eps, float dt) {
  int p_id = blockIdx.x * blockDim.x + threadIdx.x;
  if (p_id >= num_particles)
    return;

  int stack[64];
  int stack_idx = 0;

  stack[stack_idx++] = 0;

  float px = x[p_id];
  float py = y[p_id];
  float ax = 0.f, ay = 0.f;

  while (stack_idx > 0) {
    int node_idx = stack[--stack_idx];
    QuadTreeNode node = nodes[node_idx];

    if (node.total_mass == 0)
      continue;

    // leaf node
    if (node.body != -1) {
      if (p_id != node.body) {
        float dx = x[node.body] - px;
        float dy = y[node.body] - py;
        float dist_sq = dx * dx + dy * dy;
        float r = sqrt(dist_sq + eps * eps);

        float a = (G * mass[node.body]) / (r * r);
        ax += a * (dx / r);
        ay += a * (dy / r);
      }
    } else {
      float dx = node.com_x - px;
      float dy = node.com_y - py;
      float dist_sq = dx * dx + dy * dy;
      float d = sqrt(dist_sq + eps * eps);
      float width = node.boundary.halfDim * 2.0f;

      if (width / d < theta) {
        float a = (G * node.total_mass) / (d * d);
        ax += a * (dx / d);
        ay += a * (dy / d);
      } else {
        if (node.NW != -1)
          stack[stack_idx++] = node.NW;
        if (node.NE != -1)
          stack[stack_idx++] = node.NE;
        if (node.SW != -1)
          stack[stack_idx++] = node.SW;
        if (node.SE != -1)
          stack[stack_idx++] = node.SE;
      }
    }
  }

  v_x[p_id] += ax * dt;
  v_y[p_id] += ay * dt;
}
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <sys/types.h>
#include <vector>

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

struct Particles {
  int num;
  float *x;
  float *y;

  float *v_x;
  float *v_y;

  int *c_r;
  int *c_g;
  int *c_b;
  int *c_a;

  float *size;
  float *mass;

  Particles(int n) {
    num = n;
    x = new float[n];
    y = new float[n];
    v_x = new float[n];
    v_y = new float[n];
    c_r = new int[n];
    c_g = new int[n];
    c_b = new int[n];
    c_a = new int[n];
    size = new float[n];
    mass = new float[n];
  }

  vec2 pos(int id) const { return vec2(this->x[id], this->y[id]); }

  float dist(int id1, int id2) {
    return sqrt((x[id1] - x[id2]) * (x[id1] - x[id2]) +
                (y[id1] - y[id2]) * (y[id1] - y[id2]));
  }

  ~Particles() {
    delete[] x;
    delete[] y;
    delete[] v_x;
    delete[] v_y;
    delete[] c_r;
    delete[] c_g;
    delete[] c_b;
    delete[] c_a;
    delete[] size;
    delete[] mass;
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

      // True distance for collision (without softening epsilon)
      float actual_dist = std::sqrt(dist_sq);
      float min_dist = particles.size[p_id] + particles.size[root->body];

      // Collision Detection & Resolution
      if (actual_dist < min_dist && actual_dist > 0.0f) {

        // Ensure collision is only processed once per pair of particles

        float overlap = min_dist - actual_dist;
        float nx = dx / actual_dist; // Normal X
        float ny = dy / actual_dist; // Normal Y

        float m1 = particles.mass[p_id];
        float m2 = particles.mass[root->body];
        float total_mass = m1 + m2;

        // 1. Positional Correction (prevents sinking/sticking)
        float correction_x = nx * overlap;
        float correction_y = ny * overlap;

        particles.x[p_id] -= correction_x * (m2 / total_mass);
        particles.y[p_id] -= correction_y * (m2 / total_mass);

        // 2. Velocity Resolution (Elastic Bounce)
        float dvx = particles.v_x[root->body] - particles.v_x[p_id];
        float dvy = particles.v_y[root->body] - particles.v_y[p_id];

        float vel_along_normal = dvx * nx + dvy * ny;

        // Apply impulse only if they are moving towards each other
        if (vel_along_normal < 0) {
          float restitution = 0.01f; // 1.0 = perfectly bouncy, 0.0 = clay
          float j = -(1.0f + restitution) * vel_along_normal;
          j /= (1.0f / m1 + 1.0f / m2);

          float impulse_x = j * nx;
          float impulse_y = j * ny;

          particles.v_x[p_id] -= impulse_x / m1;
          particles.v_y[p_id] -= impulse_y / m1;
        }
      }
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

  ~QuadTree() {
    delete NW;
    delete NE;
    delete SW;
    delete SE;
  }
};
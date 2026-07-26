#pragma once

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
  float *c_a;

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
    c_a = new float[n];
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

  int body = -1; // stores id of the particles
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

      // adding body somewhere else
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

    // empty quadrant
    if (!root || root->total_mass == 0)
      return;

    // leaf node
    if (root->body != -1) {
      if (p_id == root->body)
        return;

      float dx = particles.x[root->body] - particles.x[p_id];
      float dy = particles.y[root->body] - particles.y[p_id];

      float r = std::sqrt(dx * dx + dy * dy + eps * eps);
      float a = (G * particles.mass[root->body]) / (r * r);

      float ax = a * (dx / r);
      float ay = a * (dy / r);

      particles.v_x[p_id] += ax * dt;
      particles.v_y[p_id] += ay * dt;
      return;
    }

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
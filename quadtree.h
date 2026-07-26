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

  constexpr static int node_capacity = 4;
  AABB boundary;

  int size = 0;
  int points[node_capacity]; // stores id of the particles

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

    if (size < node_capacity && NW == nullptr) {
      points[size++] = p_id;
      return true;
    }

    if (NW == nullptr) {
      subdivide();
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

  void queryRange(AABB Range, std::vector<int> &in_range,
                  const Particles &particles) {

    if (!boundary.intersectsAABB(Range))
      return;

    for (int i = 0; i < size; i++) {
      if (Range.containsPoint(vec2(particles.pos(points[i])))) {
        in_range.push_back(points[i]);
      }
    }

    if (!NW)
      return;

    NW->queryRange(Range, in_range, particles);
    NE->queryRange(Range, in_range, particles);
    SW->queryRange(Range, in_range, particles);
    SE->queryRange(Range, in_range, particles);

    return;
  }

  ~QuadTree() {
    delete NW;
    delete NE;
    delete SW;
    delete SE;
  }
};
#pragma once

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
  vec2 points[node_capacity];

  QuadTree *NW;
  QuadTree *NE;
  QuadTree *SW;
  QuadTree *SE;

  QuadTree(AABB boundary_): QuadTree() {
    boundary = boundary_;
  }

  QuadTree() {
    NW = nullptr;
    NE = nullptr;
    SW = nullptr;
    SE = nullptr;
  }

  bool insert(vec2 p) {

    if (!boundary.containsPoint(p))
      return false;

    if (size < node_capacity && NW == nullptr) {
      points[size++] = p;
      return true;
    }

    if (NW == nullptr) {
      subdivide();
    }

    return NW->insert(p) || NE->insert(p) || SW->insert(p) || SE->insert(p);
  }

  void subdivide() {
    vec2 center = boundary.center;
    float new_h = boundary.halfDim / 2.0f;
    NW = new QuadTree(AABB(center.x - new_h, center.y - new_h, new_h));
    NE = new QuadTree(AABB(center.x + new_h, center.y - new_h, new_h));
    SW = new QuadTree(AABB(center.x - new_h, center.y + new_h, new_h));
    SE = new QuadTree(AABB(center.x + new_h, center.y + new_h, new_h));
  }

  void queryRange(AABB Range, std::vector<vec2> &in_range) {

    if (!boundary.intersectsAABB(Range))
      return;

    for (int i = 0; i < size; i++) {
      if (Range.containsPoint(points[i])) {
        in_range.push_back(points[i]);
      }
    }

    if (!NW)
      return;

    NW->queryRange(Range, in_range);
    NE->queryRange(Range, in_range);
    SW->queryRange(Range, in_range);
    SE->queryRange(Range, in_range);

    return;
  }

  ~QuadTree() {
    delete NW;
    delete NE;
    delete SW;
    delete SE;
  }
};
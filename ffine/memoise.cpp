#include <iostream>
#include "memoise.hpp"


std::vector<double> return_cppvec(double arg, void f(double, double*),
                                  int vec_size) {
  std::vector<double> res(vec_size, 7);
  f(arg, &res[0]);
  return res; }

double memoise_project(double arg, int k, void f(double, double*), int vec_size,
                       cache_type *cache) {
  return project(arg, k, memoise, return_cppvec, cache, f, vec_size); }

void *c_build_cache() {
  cache_type *cache;
  cache = new cache_type;
  return static_cast<void *>(cache); }

void c_destroy_cache(void *cache) {
  delete static_cast<cache_type *>(cache); }

double c_memoise_project(double arg, int k, void f(double, double*),
                         int vec_size, void *cache) {
  return memoise_project(arg, k, f, vec_size,
                         static_cast<cache_type *>(cache)); }

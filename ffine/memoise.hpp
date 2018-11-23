#include <unordered_map>
#include <vector>


extern "C" {
  void *c_build_cache();

  void c_destroy_cache(void *cache);

  double c_memoise_project(double arg, int k, void f(double, double*),
                           int vec_size, void *cache); }


typedef std::unordered_map<double, std::vector<double> > cache_type;

template <typename R, typename A, typename... Ps>
R memoise(A arg, R f(A, Ps...), std::unordered_map<A, R> *cache, Ps... params) {
  auto memoised = cache->find(arg);
  return memoised == cache->end() ? (*cache)[arg] = f(arg, params...)
                                  : memoised->second; }
template <typename A, typename... Ps>
double project(A arg, int k, std::vector<double> f(A, Ps...), Ps... params) {
  return f(arg, params...)[k]; }

double memoise_project(double arg, int k, void f(double, double*), int vec_size,
                       cache_type *cache);

#include "common_config.hpp"
#include "vector_crud.hpp"
#include "deque_crud.hpp"
#include "unordered_map_crud.hpp"

#ifndef VEC_BYIT
#define VEC_BYIT(N)                                      \
NONIUS_BENCHMARK(STRING(VEC_BYIT_##N),                   \
        [](nonius::chronometer meter) {                  \
    vector_crud<double> vec_crud(N);                     \
    vec_crud.emplace_n_element_at_a_time(N);             \
    double sum;                                          \
    meter.measure([&vec_crud, &sum]                      \
            { sum = vec_crud.traverse_by_iterator(); }); \
    std::cerr << sum;                                    \
})
#endif

#ifndef DEQ_BYIT
#define DEQ_BYIT(N)                                      \
NONIUS_BENCHMARK(STRING(DEQ_BYIT_##N),                   \
        [](nonius::chronometer meter) {                  \
    deque_crud<double> deq_crud;                         \
    deq_crud.emplace_n_element_at_a_time(N);             \
    double sum;                                          \
    meter.measure([&deq_crud, &sum]                      \
            { sum = deq_crud.traverse_by_iterator(); }); \
    std::cerr << sum;                                    \
})
#endif

#ifndef MAP_BYIT
#define MAP_BYIT(N)                                      \
NONIUS_BENCHMARK(STRING(MAP_BYIT_##N),                   \
        [](nonius::chronometer meter) {                  \
    unordered_map_crud<double> map_crud(N);              \
    map_crud.emplace_n_element_at_a_time(N);             \
    double sum;                                          \
    meter.measure([&map_crud, &sum]                      \
            { sum = map_crud.traverse_by_iterator(); }); \
    std::cerr << sum;                                    \
})
#endif

VEC_BYIT(10)
VEC_BYIT(30)
VEC_BYIT(50)
VEC_BYIT(100)
VEC_BYIT(250)
VEC_BYIT(500)
VEC_BYIT(800)
VEC_BYIT(1000)
VEC_BYIT(1500)
VEC_BYIT(3000)
VEC_BYIT(5000)
VEC_BYIT(8000)
VEC_BYIT(10000)

DEQ_BYIT(10)
DEQ_BYIT(30)
DEQ_BYIT(50)
DEQ_BYIT(100)
DEQ_BYIT(250)
DEQ_BYIT(500)
DEQ_BYIT(800)
DEQ_BYIT(1000)
DEQ_BYIT(1500)
DEQ_BYIT(3000)
DEQ_BYIT(5000)
DEQ_BYIT(8000)
DEQ_BYIT(10000)

MAP_BYIT(10)
MAP_BYIT(30)
MAP_BYIT(50)
MAP_BYIT(100)
MAP_BYIT(250)
MAP_BYIT(500)
MAP_BYIT(800)
MAP_BYIT(1000)
MAP_BYIT(1500)
MAP_BYIT(3000)
MAP_BYIT(5000)
MAP_BYIT(8000)
MAP_BYIT(10000)

#include "common_config.hpp"
#include "vector_crud.hpp"
#include "deque_crud.hpp"
#include "unordered_map_crud.hpp"

#ifndef VEC_FOR
#define VEC_FOR(N)                                   \
NONIUS_BENCHMARK(STRING(VEC_FOR_##N),                \
        [](nonius::chronometer meter) {              \
    vector_crud<double> vec_crud(N);                 \
    vec_crud.emplace_n_element_at_a_time(N);         \
    double sum;                                      \
    meter.measure([&vec_crud, &sum]                  \
            { sum = vec_crud.traverse_foreach(); }); \
    std::cerr << sum;                                \
})
#endif

#ifndef DEQ_FOR
#define DEQ_FOR(N)                                   \
NONIUS_BENCHMARK(STRING(DEQ_FOR_##N),                \
        [](nonius::chronometer meter) {              \
    deque_crud<double> deq_crud;                     \
    deq_crud.emplace_n_element_at_a_time(N);         \
    double sum;                                      \
    meter.measure([&deq_crud, &sum]                  \
            { sum = deq_crud.traverse_foreach(); }); \
    std::cerr << sum;                                \
})
#endif

#ifndef MAP_FOR
#define MAP_FOR(N)                                   \
NONIUS_BENCHMARK(STRING(MAP_FOR_##N),                \
        [](nonius::chronometer meter) {              \
    unordered_map_crud<double> map_crud(N);          \
    map_crud.emplace_n_element_at_a_time(N);         \
    double sum;                                      \
    meter.measure([&map_crud, &sum]                  \
            { sum = map_crud.traverse_foreach(); }); \
    std::cerr << sum;                                \
})
#endif

VEC_FOR(10)
VEC_FOR(30)
VEC_FOR(50)
VEC_FOR(100)
VEC_FOR(250)
VEC_FOR(500)
VEC_FOR(800)
VEC_FOR(1000)
VEC_FOR(1500)
VEC_FOR(3000)
VEC_FOR(5000)
VEC_FOR(8000)
VEC_FOR(10000)

DEQ_FOR(10)
DEQ_FOR(30)
DEQ_FOR(50)
DEQ_FOR(100)
DEQ_FOR(250)
DEQ_FOR(500)
DEQ_FOR(800)
DEQ_FOR(1000)
DEQ_FOR(1500)
DEQ_FOR(3000)
DEQ_FOR(5000)
DEQ_FOR(8000)
DEQ_FOR(10000)

MAP_FOR(10)
MAP_FOR(30)
MAP_FOR(50)
MAP_FOR(100)
MAP_FOR(250)
MAP_FOR(500)
MAP_FOR(800)
MAP_FOR(1000)
MAP_FOR(1500)
MAP_FOR(3000)
MAP_FOR(5000)
MAP_FOR(8000)
MAP_FOR(10000)

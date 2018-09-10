#include "common_config.hpp"
#include "vector_crud.hpp"
#include "deque_crud.hpp"
#include "unordered_map_crud.hpp"

#ifndef VEC_BYID
#define VEC_BYID(N)                                \
NONIUS_BENCHMARK(STRING(VEC_BYID_##N),             \
        [](nonius::chronometer meter) {            \
    vector_crud<double> vec_crud(N);               \
    vec_crud.emplace_n_element_at_a_time(N);       \
    double sum;                                    \
    meter.measure([&vec_crud, &sum]                \
            { sum = vec_crud.traverse_by_id(); }); \
    std::cerr << sum;                              \
})
#endif

#ifndef DEQ_BYID
#define DEQ_BYID(N)                                \
NONIUS_BENCHMARK(STRING(DEQ_BYID_##N),             \
        [](nonius::chronometer meter) {            \
    deque_crud<double> deq_crud;                   \
    deq_crud.emplace_n_element_at_a_time(N);       \
    double sum;                                    \
    meter.measure([&deq_crud, &sum]                \
            { sum = deq_crud.traverse_by_id(); }); \
    std::cerr << sum;                              \
})
#endif

#ifndef MAP_BYID
#define MAP_BYID(N)                                \
NONIUS_BENCHMARK(STRING(MAP_BYID_##N),             \
        [](nonius::chronometer meter) {            \
    unordered_map_crud<double> map_crud(N);        \
    map_crud.emplace_n_element_at_a_time(N);       \
    double sum;                                    \
    meter.measure([&map_crud, &sum]                \
            { sum = map_crud.traverse_by_id(); }); \
    std::cerr << sum;                              \
})
#endif

VEC_BYID(10)
VEC_BYID(30)
VEC_BYID(50)
VEC_BYID(100)
VEC_BYID(250)
VEC_BYID(500)
VEC_BYID(800)
VEC_BYID(1000)
VEC_BYID(1500)
VEC_BYID(3000)
VEC_BYID(5000)
VEC_BYID(8000)
VEC_BYID(10000)

DEQ_BYID(10)
DEQ_BYID(30)
DEQ_BYID(50)
DEQ_BYID(100)
DEQ_BYID(250)
DEQ_BYID(500)
DEQ_BYID(800)
DEQ_BYID(1000)
DEQ_BYID(1500)
DEQ_BYID(3000)
DEQ_BYID(5000)
DEQ_BYID(8000)
DEQ_BYID(10000)

MAP_BYID(10)
MAP_BYID(30)
MAP_BYID(50)
MAP_BYID(100)
MAP_BYID(250)
MAP_BYID(500)
MAP_BYID(800)
MAP_BYID(1000)
MAP_BYID(1500)
MAP_BYID(3000)
MAP_BYID(5000)
MAP_BYID(8000)
MAP_BYID(10000)

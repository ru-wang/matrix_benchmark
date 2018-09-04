#pragma once

#include "./common_config.hpp"

template<typename element_t>
struct unordered_map_crud {
    static_assert(std::is_arithmetic<element_t>::value,
            "currently we only test arithmetic types");
    static_assert(not std::is_reference<element_t>::value,
            "reference types are not allowed");
    using element_sum_type = std::remove_const_t<element_t>;

    template<int remain, typename ele_t, typename ...args_t>
        struct emplace_recurse {
            static void impl(std::unordered_map<int, ele_t>& container,
                    ele_t& ele, args_t... args) {
                container.emplace((int)container.size(), ele);
                emplace_recurse<sizeof...(args_t)-1, args_t...>
                    ::impl(container, args...);
            }
        };
    template<typename ele_t>
        struct emplace_recurse<0, ele_t> {
            static void impl(std::unordered_map<int, ele_t>& container, ele_t& ele)
            { container.emplace((int)container.size(), ele); }
        };

    unordered_map_crud() = default;
    unordered_map_crud(int n) { _container_.reserve(n); }

    /* CRUD: create */
    template<typename ...args_t>
        void emplace(args_t... args) {
            emplace_recurse<sizeof...(args)-1, args_t...>
                ::impl(_container_, args...);
        }
    void emplace_n_element_at_a_time(int n) {
        for (int i = 0; i < n; ++i)
            _container_.emplace(i, (element_t)i);
    }

    /* CRUD: read */
    element_sum_type traverse_foreach() const {
        element_sum_type sum = 0;
        for (auto& ele : _container_)
            sum += ele.second;
        return sum;
    }
    element_sum_type traverse_by_iterator() const {
        element_sum_type sum = 0;
        for (auto it = _container_.begin(); it != _container_.end(); ++it)
            sum += it->second;
        return sum;
    }
    element_sum_type traverse_by_id() const {
        element_sum_type sum = 0;
        for (int i = 0; i < (int)_container_.size(); ++i)
            sum += _container_.at(i);
        return sum;
    }
    element_t& random_access(int id)
    { return _container_.at(id); }

    /* CRUD: update */

    /* CRUD: delete */

    std::unordered_map<int, element_t> _container_;
};

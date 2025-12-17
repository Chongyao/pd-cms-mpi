#include "timer.hpp"
#include <iostream>
#include <mutex>
namespace zcy {

TIMING& TIMING::instance()
{
    static TIMING m_instance;
    return m_instance;
}

void TIMING::__TIME_BEGIN__(const std::string& key, size_t thread_id)
{
    auto& timer_list = get_timer_list(thread_id);
    auto it = timer_list.find(key);
    if (it != timer_list.end()) {
        std::cerr << key << " Timer conflict.\n";
        // print timer_list
        for (const auto& [k, v] : timer_list) {
            std::cerr << k << " : " << v.time_since_epoch().count() << "\n";
        }
        throw -1;
    }
    timer_list.insert({ key, std::chrono::steady_clock::now() });
    // timer_list.push_back(std::chrono::steady_clock::now());
}

std::unordered_map<std::string, std::chrono::steady_clock::time_point>& TIMING::get_timer_list(size_t thread_id)
{
    if (timer_list_threads.size() <= thread_id) {
        std::cout << "resize timer list \n";
        timer_list_threads.resize(thread_id + 1);
    }
    return timer_list_threads[thread_id];
}
std::unordered_map<std::string, double>& TIMING::get_timer(size_t thread_id)
{
    if (timer_per_thread.size() <= thread_id) {
        timer_per_thread.resize(thread_id + 1);
    }
    return timer_per_thread[thread_id];
}
double TIMING::__TIME_END__(const std::string& key, size_t thread_id)
{
    auto& timer_list = get_timer_list(thread_id);
    // if (timer_list.empty()) {
    //     std::cout << "missing a __TIME__BEGIN__;" << std::endl;
    //     throw -1;
    // }
    auto it = timer_list.find(key);
    if (it == timer_list.end()) {
        std::cerr << key << " Missing a TIMING::TIC() / __TIME_BEGIN__.\n";
        throw -1;
    }
    auto start = it->second;
    timer_list.erase(it);
    auto end = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    double res = double(duration.count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den;

    return res;
}

std::unordered_map<std::string, std::chrono::steady_clock::time_point>& EXTRA_TIMING::get_timer_list()
{
    return timer_list;
}
void EXTRA_TIMING::__TIME_BEGIN__(const std::string& key)
{
    auto& timer_list = get_timer_list();
    auto it = timer_list.find(key);
    if (it != timer_list.end()) {
        std::cerr << key << " Timer conflict.\n";
        // print timer_list
        for (const auto& [k, v] : timer_list) {
            std::cerr << k << " : " << v.time_since_epoch().count() << "\n";
        }
        throw -1;
    }
    timer_list.insert({ key, std::chrono::steady_clock::now() });
    // timer_list.push_back(std::chrono::steady_clock::now());
}
double EXTRA_TIMING::__TIME_END__(const std::string& info, bool if_print)
{
    auto& timer_list = get_timer_list();
    // if (timer_list.empty()) {
    //     std::cout << "missing a __TIME__BEGIN__;" << std::endl;
    //     throw -1;
    // }
    auto it = timer_list.find(info);
    if (it == timer_list.end()) {
        std::cerr << info << " Missing a TIMING::TIC() / __TIME_BEGIN__.\n";
        throw -1;
    }
    auto start = it->second;
    // timer_list.erase(it);
    auto end = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    double res = double(duration.count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den;
    if (if_print) {
        std::cout << info << " cost " << res << "seoncds\n";
    }
    return res;
}
EXTRA_TIMING& EXTRA_TIMING::instance()
{
    static EXTRA_TIMING m_instance;
    return m_instance;
}

}

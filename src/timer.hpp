#pragma once

#include <unordered_map>
// #include <spdlog/fmt/ranges.h>
#include <fmt/ranges.h>
#include <spdlog/spdlog.h>
#include <string>

namespace zcy {
struct ScopedTimer {
    std::chrono::steady_clock::time_point start;
    double& result_in_seconds;

    explicit ScopedTimer(double& result)
        : result_in_seconds(result)
    {
        start = std::chrono::steady_clock::now();
    }

    ~ScopedTimer()
    {
        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        result_in_seconds = elapsed.count();
    }

    ScopedTimer(const ScopedTimer&) = delete;
    ScopedTimer& operator=(const ScopedTimer&) = delete;
    ScopedTimer(ScopedTimer&&) = delete;
    ScopedTimer& operator=(ScopedTimer&&) = delete;
};

struct TIMING {
private:
    std::vector<std::unordered_map<std::string, double>> timer_per_thread;
    // std::unordered_map<std::string, double> timer;
    // std::unordered_map<std::string, std::chrono::steady_clock::time_point> timer_list;
    std::vector<std::unordered_map<std::string, std::chrono::steady_clock::time_point>> timer_list_threads;

public:
    static void set_num_threads(size_t num)
    {
        instance().timer_per_thread.resize(num);
        instance().timer_list_threads.resize(num);
    }
    bool is_enable { true };
    static void ENABLE()
    {
        instance().is_enable = true;
    }
    static void DISABLE()
    {
        instance().is_enable = false;
    }

    static void TIC(const std::string& key, bool print_now = false, size_t thread_id = 0)
    {

        if (instance().is_enable) {
            instance().__TIME_BEGIN__(key, thread_id);
            if (print_now) {
                spdlog::info("TIC: {}", key);
            }
        }
    }

    static double TOC(const std::string& key, bool print_now = false, size_t thread_id = 0)
    {
        if (instance().is_enable) {
            double ret = instance().__TIME_END__(key, thread_id);
            if (print_now) {
                spdlog::info("TOC: {} cost {} seconds.", key, ret);
            }
            auto& timer = instance().get_timer(thread_id);
            if (timer.find(key) != timer.end()) {
                timer[key] += ret;
            } else {
                timer[key] = ret;
            }
            return ret;
        } else {
            return 0;
        }
    }
    static void reset(size_t thread_id = 0)
    {
        instance().get_timer(thread_id).clear();
    }
    static void print(size_t threa_id = 0)
    {
        spdlog::info("TIMING record: {}\n", instance().get_timer(threa_id));
    }
    static std::unordered_map<std::string, double> sum_all_threads()
    {
        std::unordered_map<std::string, double> res;
        for (const auto& t : instance().timer_per_thread) {
            for (const auto& [k, v] : t) {
                if (res.find(k) != res.end()) {
                    res[k] += v;
                } else {
                    res[k] = v;
                }
            }
        }
        return res;
    }
    static double get_record(const std::string& key, size_t thread_id = 0)
    {
        const auto& timer = instance().get_timer(thread_id);
        const auto it = timer.find(key);
        double t = 0;
        if (it != timer.end()) {
            t = it->second;
        }
        return t;
    }

private:
    TIMING() = default;
    static TIMING& instance();
    std::unordered_map<std::string, std::chrono::steady_clock::time_point>& get_timer_list(size_t thread_id);
    std::unordered_map<std::string, double>& get_timer(size_t thread_id);
    void __TIME_BEGIN__(const std::string& key, size_t thread_id);
    double __TIME_END__(const std::string& key, size_t thread_id);
};

struct EXTRA_TIMING {
private:
    std::unordered_map<std::string, std::chrono::steady_clock::time_point> timer_list;

public:
    static void TIC(const std::string& key)
    {
        instance().__TIME_BEGIN__(key);
    }
    static double TOC(const std::string& key, bool print_now = false)
    {
        double ret = instance().__TIME_END__(key, print_now);
        return ret;
    }

private:
    EXTRA_TIMING() = default;
    static EXTRA_TIMING& instance();
    std::unordered_map<std::string, std::chrono::steady_clock::time_point>& get_timer_list();
    void __TIME_BEGIN__(const std::string& key);
    double __TIME_END__(const std::string& info, bool if_print = true);
};

}

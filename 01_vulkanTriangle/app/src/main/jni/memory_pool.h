#pragma once
//#include "taichi/common/core.h"
#include "unified_allocator.h"
#define TI_RUNTIME_HOST
#include "mem_request.h"

#include <mutex>
#include <vector>
#include <memory>
#include <thread>

//class Program;

// A memory pool that runs on the host

class MemoryPool {
public:
    std::vector<std::unique_ptr<UnifiedAllocator>> allocators;
    static constexpr std::size_t default_allocator_size =
            1 << 30;  // 1 GB per allocator
    bool terminating, killed;
    std::mutex mut;
    std::mutex mut_allocators;
    std::unique_ptr<std::thread> th;
    int processed_tail;
    bool use_unified_memory;
//  Program *prog;

    MemRequestQueue *queue;
    void *cuda_stream{nullptr};

    MemoryPool();

    template <typename T>
    T fetch(volatile void *ptr);

    template <typename T>
    void push(volatile T *dest, const T &val);

    void *allocate(std::size_t size, std::size_t alignment);

    void set_queue(MemRequestQueue *queue);

    void daemon();

    void terminate();

    ~MemoryPool();

private:
    static constexpr bool use_cuda_stream = false;
};

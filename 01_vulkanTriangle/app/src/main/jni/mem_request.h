#pragma once

// Use relative path here for runtime compilation
#include "constants.h"


struct MemRequest {
    std::size_t size;
    std::size_t alignment;
    uint8 *ptr;
    std::size_t __padding;
};

static_assert((sizeof(MemRequest) & (sizeof(MemRequest) - 1)) == 0);

struct MemRequestQueue {
    MemRequest requests[taichi_max_num_mem_requests];
    int tail;
    int processed;
};


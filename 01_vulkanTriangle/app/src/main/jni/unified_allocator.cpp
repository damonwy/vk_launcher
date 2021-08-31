// Virtual memory allocator for CPU/GPU

#include "unified_allocator.h"
#include "virtual_memory.h"
#include <string>
#include <cassert>

UnifiedAllocator::UnifiedAllocator(std::size_t size)
        : size(size) {

//    TI_TRACE("Allocating virtual address space of size {} MB",
//             size / 1024 / 1024);
    cpu_vm = std::make_unique<VirtualMemoryAllocator>(size);
    data = (uint8 *)cpu_vm->ptr;
    assert(data != nullptr);
    assert(uint64(data) % 4096 == 0);
//  TI_ASSERT(data != nullptr);
//  TI_ASSERT(uint64(data) % 4096 == 0);

    head = data;
    tail = head + size;
//  TI_TRACE("Memory allocated. Allocation time = {:.3} s", Time::get_time() - t);
}

UnifiedAllocator::~UnifiedAllocator() {
    if (!initialized()) {
        return;
    }
}

void UnifiedAllocator::memset(unsigned char val) {
    std::memset(data, val, size);
}


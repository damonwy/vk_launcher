#pragma once

#include <vulkan_wrapper.h>
#include <stdexcept>

#pragma message("BAIL_ON_VK_BAD_RESULT uses exception")

#define BAIL_ON_VK_BAD_RESULT(result, msg) \
  do {                                     \
    if ((result) != VK_SUCCESS) {          \
      throw std::runtime_error((msg));     \
    };                                     \
  } while (0)

inline constexpr VkAllocationCallbacks *kNoVkAllocCallbacks = nullptr;


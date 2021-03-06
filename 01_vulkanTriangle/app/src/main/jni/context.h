#pragma once
#include <cstdio>
#include <cinttypes>
// Use relative path here for runtime compilation
//#include "constants.h"

using int32 = int32_t;
using uint64 = uint64_t;

using ContextArgType = long long;
struct LLVMRuntime;
constexpr int taichi_max_num_indices = 8;
constexpr int taichi_max_num_args = 8;

//template <typename T, typename G>
//T taichi_union_cast_with_different_sizes(G g) {
//  union {
//    T t;
//    G g;
//  } u;
//  u.g = g;
//  return u.t;
//}

// "Context" holds necessary data for kernel body execution, such as a pointer
// to the LLVMRuntime struct, kernel arguments, and the thread id (if on CPU).
//struct Context {
//  LLVMRuntime *runtime;
//  uint64 args[taichi_max_num_args];
//  int32 extra_args[taichi_max_num_args][taichi_max_num_indices];
//  int32 cpu_thread_id;
//
//  static constexpr size_t extra_args_size = sizeof(extra_args);
//
//  template <typename T>
//  T get_arg(int i) {
//    return taichi_union_cast_with_different_sizes<T>(args[i]);
//  }
//
//  uint64 get_arg_as_uint64(int i) {
//    return args[i];
//  }
//
//  template <typename T>
//  void set_arg(int i, T v) {
//    args[i] = taichi_union_cast_with_different_sizes<uint64>(v);
//  }
//};


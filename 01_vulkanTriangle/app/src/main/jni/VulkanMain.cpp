//
// Created by Daosheng Mu on 8/8/20.
//

#include "VulkanMain.h"
#include <android/log.h>
#include <android_native_app_glue.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>
#include <stdexcept>
#include <cstdlib>
#include <spdlog/fmt/fmt.h>
#include "cassert"
#include "vulkan_wrapper.h"
#include "VulkanRenderer.h"
#include "Logger.h"

#include "vulkan_api.h"
#include "kernel_utils.h"
#include "common.h"
#include "runtime.h"
#include "vulkan_simple_memory_pool.h"
#include "constants.h"
#include "memory_pool.h"

template <typename T, typename G>
T taichi_union_cast_with_different_sizes(G g) {
  union {
      T t;
      G g;
  } u;
  u.g = g;
  return u.t;
}

struct Context {
//  LLVMRuntime *runtime;
    uint64 args[taichi_max_num_args];
    int32 extra_args[taichi_max_num_args][taichi_max_num_indices];
//    int32 cpu_thread_id;

    static constexpr size_t extra_args_size = sizeof(extra_args);

    template <typename T>
    T get_arg(int i) {
      return taichi_union_cast_with_different_sizes<T>(args[i]);
    }

    uint64 get_arg_as_uint64(int i) {
      return args[i];
    }

    template <typename T>
    void set_arg(int i, T v) {
      args[i] = taichi_union_cast_with_different_sizes<uint64>(v);
    }
};


//class Kernel {
//public:
//    std::string name;
//    struct Arg {
//        DataType dt;
//        bool is_external_array;
//        std::size_t size;
//
//        explicit Arg(const DataType &dt = PrimitiveType::unknown,
//                     bool is_external_array = false,
//                     std::size_t size = 0)
//                : dt(dt), is_external_array(is_external_array), size(size) {
//        }
//    };
//
//    struct Ret {
//        DataType dt;
//
//        explicit Ret(const DataType &dt = PrimitiveType::unknown) : dt(dt) {
//        }
//    };
//
//    std::vector<Arg> args;
//    std::vector<Ret> rets;
//    // TODO: Give "Context" a more specific name.
//    class LaunchContextBuilder {
//    public:
//        LaunchContextBuilder(Kernel *kernel, Context *ctx);
//        explicit LaunchContextBuilder(Kernel *kernel);
//
//        LaunchContextBuilder(LaunchContextBuilder &&) = default;
//        LaunchContextBuilder &operator=(LaunchContextBuilder &&) = default;
//        LaunchContextBuilder(const LaunchContextBuilder &) = delete;
//        LaunchContextBuilder &operator=(const LaunchContextBuilder &) = delete;
//
//        void set_arg_float(int arg_id, float64 d);
//
//        void set_arg_int(int arg_id, int64 d);
//
//        void set_extra_arg_int(int i, int j, int32 d);
//
//        void set_arg_external_array(int arg_id, uint64 ptr, uint64 size);
//
//        // Sets the |arg_id|-th arg in the context to the bits stored in |d|.
//        // This ignores the underlying kernel's |arg_id|-th arg type.
//        void set_arg_raw(int arg_id, uint64 d);
//
//        Context &get_context();
//
//    private:
//        Kernel *kernel_;
//        std::unique_ptr<Context> owned_ctx_;
//        // |ctx_| *almost* always points to |owned_ctx_|. However, it is possible
//        // that the caller passes a Context pointer externally. In that case,
//        // |owned_ctx_| will be nullptr.
//        // Invariant: |ctx_| will never be nullptr.
//        Context *ctx_;
//    };
//
//};
//
//Kernel::LaunchContextBuilder::LaunchContextBuilder(Kernel *kernel, Context *ctx)
//        : kernel_(kernel), owned_ctx_(nullptr), ctx_(ctx) {
//}
//
//Kernel::LaunchContextBuilder::LaunchContextBuilder(Kernel *kernel)
//        : kernel_(kernel),
//          owned_ctx_(std::make_unique<Context>()),
//          ctx_(owned_ctx_.get()) {
//}
//
//void Kernel::LaunchContextBuilder::set_arg_float(int arg_id, float64 d) {
////  TI_ASSERT_INFO(!kernel_->args[arg_id].is_external_array,
////                 "Assigning scalar value to external(numpy) array argument is "
////                 "not allowed.");
//
////  ActionRecorder::get_instance().record(
////      "set_kernel_arg_float64",
////      {ActionArg("kernel_name", kernel_->name), ActionArg("arg_id", arg_id),
////       ActionArg("val", d)});
//
//  auto dt = kernel_->args[arg_id].dt;
//  if (dt->is_primitive(PrimitiveTypeID::f32)) {
//    ctx_->set_arg(arg_id, (float32)d);
//  } else if (dt->is_primitive(PrimitiveTypeID::f64)) {
//    ctx_->set_arg(arg_id, (float64)d);
//  } else if (dt->is_primitive(PrimitiveTypeID::i32)) {
//    ctx_->set_arg(arg_id, (int32)d);
//  } else if (dt->is_primitive(PrimitiveTypeID::i64)) {
//    ctx_->set_arg(arg_id, (int64)d);
//  } else if (dt->is_primitive(PrimitiveTypeID::i8)) {
//    ctx_->set_arg(arg_id, (int8)d);
//  } else if (dt->is_primitive(PrimitiveTypeID::i16)) {
//    ctx_->set_arg(arg_id, (int16)d);
//  } else if (dt->is_primitive(PrimitiveTypeID::u8)) {
//    ctx_->set_arg(arg_id, (uint8)d);
//  } else if (dt->is_primitive(PrimitiveTypeID::u16)) {
//    ctx_->set_arg(arg_id, (uint16)d);
//  } else if (dt->is_primitive(PrimitiveTypeID::u32)) {
//    ctx_->set_arg(arg_id, (uint32)d);
//  } else if (dt->is_primitive(PrimitiveTypeID::u64)) {
//    ctx_->set_arg(arg_id, (uint64)d);
//  } else {
////    TI_NOT_IMPLEMENTED
//  }
//}
//
//void Kernel::LaunchContextBuilder::set_arg_int(int arg_id, int64 d) {
////  TI_ASSERT_INFO(!kernel_->args[arg_id].is_external_array,
////                 "Assigning scalar value to external(numpy) array argument is "
////                 "not allowed.");
//
////  ActionRecorder::get_instance().record(
////      "set_kernel_arg_int64",
////      {ActionArg("kernel_name", kernel_->name), ActionArg("arg_id", arg_id),
////       ActionArg("val", d)});
//
//  auto dt = kernel_->args[arg_id].dt;
//  if (dt->is_primitive(PrimitiveTypeID::i32)) {
//    ctx_->set_arg(arg_id, (int32)d);
//  } else if (dt->is_primitive(PrimitiveTypeID::i64)) {
//    ctx_->set_arg(arg_id, (int64)d);
//  } else if (dt->is_primitive(PrimitiveTypeID::i8)) {
//    ctx_->set_arg(arg_id, (int8)d);
//  } else if (dt->is_primitive(PrimitiveTypeID::i16)) {
//    ctx_->set_arg(arg_id, (int16)d);
//  } else if (dt->is_primitive(PrimitiveTypeID::u8)) {
//    ctx_->set_arg(arg_id, (uint8)d);
//  } else if (dt->is_primitive(PrimitiveTypeID::u16)) {
//    ctx_->set_arg(arg_id, (uint16)d);
//  } else if (dt->is_primitive(PrimitiveTypeID::u32)) {
//    ctx_->set_arg(arg_id, (uint32)d);
//  } else if (dt->is_primitive(PrimitiveTypeID::u64)) {
//    ctx_->set_arg(arg_id, (uint64)d);
//  } else if (dt->is_primitive(PrimitiveTypeID::f32)) {
//    ctx_->set_arg(arg_id, (float32)d);
//  } else if (dt->is_primitive(PrimitiveTypeID::f64)) {
//    ctx_->set_arg(arg_id, (float64)d);
//  } else {
////    TI_INFO(dt->to_string());
////    TI_NOT_IMPLEMENTED
//  }
//}
//
//void Kernel::LaunchContextBuilder::set_extra_arg_int(int i, int j, int32 d) {
//  ctx_->extra_args[i][j] = d;
//}
//
//void Kernel::LaunchContextBuilder::set_arg_external_array(int arg_id,
//                                                          uint64 ptr,
//                                                          uint64 size) {
////  TI_ASSERT_INFO(
////      kernel_->args[arg_id].is_external_array,
////      "Assigning external(numpy) array to scalar argument is not allowed.");
////
////  ActionRecorder::get_instance().record(
////      "set_kernel_arg_ext_ptr",
////      {ActionArg("kernel_name", kernel_->name), ActionArg("arg_id", arg_id),
////       ActionArg("address", fmt::format("0x{:x}", ptr)),
////       ActionArg("array_size_in_bytes", (int64)size)});
//
//  kernel_->args[arg_id].size = size;
//  ctx_->set_arg(arg_id, ptr);
//}
//
//void Kernel::LaunchContextBuilder::set_arg_raw(int arg_id, uint64 d) {
////  TI_ASSERT_INFO(!kernel_->args[arg_id].is_external_array,
////                 "Assigning scalar value to external(numpy) array argument is "
////                 "not allowed.");
////
////  if (!kernel_->is_evaluator) {
////    ActionRecorder::get_instance().record(
////        "set_arg_raw",
////        {ActionArg("kernel_name", kernel_->name), ActionArg("arg_id", arg_id),
////         ActionArg("val", (int64)d)});
////  }
//  ctx_->set_arg<uint64>(arg_id, d);
//}

//Context &Kernel::LaunchContextBuilder::get_context() {
////  ctx_->runtime = static_cast<LLVMRuntime *>(kernel_->program->llvm_runtime);
//  return *ctx_;
//}

class HostDeviceContextBlitter {
public:
    HostDeviceContextBlitter(const KernelContextAttributes *ctx_attribs,
                             Context *host_ctx,
                             uint64_t *host_result_buffer,
                             VkBufferWithMemory *device_buffer)
            : ctx_attribs_(ctx_attribs),
              host_ctx_(host_ctx),
              host_result_buffer_(host_result_buffer),
              device_buffer_(device_buffer) {
    }

    void host_to_device() {
      if (ctx_attribs_->empty()) {
        return;
      }
      auto mapped = device_buffer_->map_mem();
      char *const device_base = reinterpret_cast<char *>(mapped.data());

#define TO_DEVICE(short_type, type)                         \
  else if (dt->is_primitive(PrimitiveTypeID::short_type)) { \
    auto d = host_ctx_->get_arg<type>(i);                   \
    std::memcpy(device_ptr, &d, sizeof(d));                 \
  }
      for (int i = 0; i < ctx_attribs_->args().size(); ++i) {
        const auto &arg = ctx_attribs_->args()[i];
//      const auto dt = arg.dt;
        const auto dt = PrimitiveType::f32;

        char *device_ptr = device_base + arg.offset_in_mem;
        if (arg.is_array) {
          const void *host_ptr = host_ctx_->get_arg<void *>(i);
          std::memcpy(device_ptr, host_ptr, arg.stride);
        }
//      TO_DEVICE(i32, int32)
        if (dt->is_primitive(PrimitiveTypeID::i32)) {
          auto d = host_ctx_->get_arg<int32>(i);
          std::memcpy(device_ptr, &d, sizeof(d));

        }
//      TO_DEVICE(u32, uint32)
        if (dt->is_primitive(PrimitiveTypeID::u32)) {
          auto d = host_ctx_->get_arg<uint32>(i);
          std::memcpy(device_ptr, &d, sizeof(d));

        }
//      TO_DEVICE(f32, float32)
        if (dt->is_primitive(PrimitiveTypeID::f32)) {
          auto d = host_ctx_->get_arg<float32>(i);
          std::memcpy(device_ptr, &d, sizeof(d));

        }
        else {
//        TI_ERROR("Vulkan does not support arg type={}", data_type_name(arg.dt));
        }
      }
      char *device_ptr = device_base + ctx_attribs_->extra_args_mem_offset();
      std::memcpy(device_ptr, host_ctx_->extra_args,
                  ctx_attribs_->extra_args_bytes());
#undef TO_DEVICE
    }

    void device_to_host() {
      if (ctx_attribs_->empty()) {
        return;
      }
      auto mapped = device_buffer_->map_mem();
      char *const device_base = reinterpret_cast<char *>(mapped.data());
#define TO_HOST(short_type, type)                           \
  else if (dt->is_primitive(PrimitiveTypeID::short_type)) { \
    const type d = *reinterpret_cast<type *>(device_ptr);   \
    host_result_buffer_[i] =                                \
        taichi_union_cast_with_different_sizes<uint64>(d);  \
  }
      for (int i = 0; i < ctx_attribs_->args().size(); ++i) {
        const auto &arg = ctx_attribs_->args()[i];
        char *device_ptr = device_base + arg.offset_in_mem;
        if (arg.is_array) {
          void *host_ptr = host_ctx_->get_arg<void *>(i);
          std::memcpy(host_ptr, device_ptr, arg.stride);
        }
      }
      for (int i = 0; i < ctx_attribs_->rets().size(); ++i) {
        // Note that we are copying the i-th return value on Metal to the i-th
        // *arg* on the host context.
        const auto &ret = ctx_attribs_->rets()[i];
        char *device_ptr = device_base + ret.offset_in_mem;
//      const auto dt = ret.dt;
        const auto dt = PrimitiveType::f32;

        if (ret.is_array) {
          void *host_ptr = host_ctx_->get_arg<void *>(i);
          std::memcpy(host_ptr, device_ptr, ret.stride);
        }
        TO_HOST(i32, int32)
        TO_HOST(u32, uint32)
        TO_HOST(f32, float32)
        else {
//        TI_ERROR("Vulkan does not support return value type={}",
//                 data_type_name(ret.dt));
        }
      }
#undef TO_HOST
    }

    static std::unique_ptr<HostDeviceContextBlitter> maybe_make(
            const KernelContextAttributes *ctx_attribs,
            Context *host_ctx,
            uint64_t *host_result_buffer,
            VkBufferWithMemory *device_buffer) {
      if (ctx_attribs->empty()) {
        return nullptr;
      }
      return std::make_unique<HostDeviceContextBlitter>(ctx_attribs, host_ctx, host_result_buffer, device_buffer);
    }

private:
    const KernelContextAttributes *const ctx_attribs_;
    Context *const host_ctx_;
    uint64_t *const host_result_buffer_;
    VkBufferWithMemory *const device_buffer_;
};



class KernelHandle {
public:
//  friend class Impl;
    int id_ = -1;
};

struct Params {
    // CompiledSNodeStructs compiled_snode_structs;
    uint64_t *host_result_buffer = nullptr;
    // int root_id;
//    const SNodeDescriptorsMap *snode_descriptors = nullptr;
};

class ClearBufferCommandBuilder : private VulkanCommandBuilder {
public:
    using VulkanCommandBuilder::VulkanCommandBuilder;

    VkCommandBuffer build(const std::vector<VkBuffer> &buffers) {
      for (auto b : buffers) {
        vkCmdFillBuffer(command_buffer_, b, /*dstOffset=*/0,
                /*size=*/VK_WHOLE_SIZE,
                /*data=*/0);
      }
      return VulkanCommandBuilder::build();
    }
};

class Test {
public:
    Test(const Params &params, android_app* app):
            host_result_buffer_(params.host_result_buffer){
      ManagedVulkanDevice::Params mvd_params;
      mvd_params.api_version = VulkanEnvSettings::kApiVersion();
      managed_device_ = std::make_unique<ManagedVulkanDevice>(mvd_params, app);
      stream_ = std::make_unique<VulkanStream>(managed_device_->device());

      init_memory_pool(params);
      init_vk_buffers();
    }

    void run() {
//        auto kernel_run = std::make_unique<Kernel>();
      auto host_ctx = new Context();
      host_ctx->set_arg(0, (float32)9.8);

      load_kernels();
      std::cout << "Load kernels successfully!... " << std::endl;
      std::vector<KernelHandle> vec_handles;
      for(auto const &rp: register_params_){
        auto handle = register_taichi_kernel(rp);
        vec_handles.push_back(handle);
        launch_kernel(handle, host_ctx);
      }
      auto map_root = root_buffer_->map_mem();
      char *const root_base = reinterpret_cast<char *>(map_root.data());

//        const unsigned char *charBuffer = (unsigned char *) root_base;
      //Then we create the vector (named vectorBuffer) by copying the contents of charBuffer to the vector
      std::vector<unsigned char> vectorBuffer(root_base, root_base + sizeof(uint32_t)*8);

      char res[32];
      std::memcpy(res, root_base, sizeof(uint32_t)*8);

      std::vector<float> result(8);
      // Copy data from the char-vector to a new uint32_t-vector
      memcpy(result.data(), vectorBuffer.data(), 32);
      for (int i = 0; i < 8; ++i){
        std::cout << result[i] << std::endl;

      }
      for(auto const &handle: vec_handles){
        free_command_buffer(handle);
      }

    }

    ~Test(){
      global_tmps_buffer_.reset();
      root_buffer_.reset();
      dev_local_memory_pool_.reset();
    }

    void load_kernels(){
      RegisterParams reg_params0;
      reg_params0.kernel_attribs.name = "init_c6_0_k0004_vk";
      reg_params0.kernel_attribs.is_jit_evaluator = false;
      TaskAttributes t00;
      t00.advisory_total_num_threads = 1;
      t00.advisory_num_threads_per_group = 1;
      t00.task_type = 0;
      t00.buffer_binds.push_back({BufferEnum::Root,0});
      t00.buffer_binds.push_back({BufferEnum::GlobalTmps,1});
//        t0.buffer_binds.push_back({BufferEnum::Context,2});
      t00.name = "init_c6_0_k0004_vk_t00";

      TaskAttributes t01;
      t01.advisory_total_num_threads = 5;
      t01.advisory_num_threads_per_group = 128;
      t01.task_type = 1;
      t01.buffer_binds.push_back({BufferEnum::Root,0});
      t01.buffer_binds.push_back({BufferEnum::GlobalTmps,1});
//        t0.buffer_binds.push_back({BufferEnum::Context,2});
      t01.name = "init_c6_0_k0004_vk_t01";
      t01.range_for_attribs = std::make_optional<TaskAttributes::RangeForAttributes>();
      t01.range_for_attribs->begin = 0;
      t01.range_for_attribs->end = 5;
      t01.range_for_attribs->const_begin = true;
      t01.range_for_attribs->const_end = true;

      reg_params0.kernel_attribs.tasks_attribs.push_back(t00);
      reg_params0.kernel_attribs.tasks_attribs.push_back(t01);

      reg_params0.kernel_attribs.ctx_attribs.args_bytes_=0;
      reg_params0.kernel_attribs.ctx_attribs.rets_bytes_=0;
      reg_params0.kernel_attribs.ctx_attribs.extra_args_bytes_=256;
      register_params_.push_back(std::move(reg_params0));

      RegisterParams reg_params1;
      reg_params1.kernel_attribs.name = "substep_c4_0_k0006_vk";
      reg_params1.kernel_attribs.is_jit_evaluator = false;
      TaskAttributes t10;
      t10.advisory_total_num_threads = 256;
      t10.advisory_num_threads_per_group = 128;
      t10.task_type = 1;
      t10.buffer_binds.push_back({BufferEnum::Root,0});
      t10.buffer_binds.push_back({BufferEnum::GlobalTmps,1});
      t10.buffer_binds.push_back({BufferEnum::Context,2});
      t10.name = "substep_c4_0_k0006_vk_t00";
      t10.range_for_attribs = std::make_optional<TaskAttributes::RangeForAttributes>();
      t10.range_for_attribs->begin = 0;
      t10.range_for_attribs->end = 256;
      t10.range_for_attribs->const_begin = true;
      t10.range_for_attribs->const_end = true;
      reg_params1.kernel_attribs.tasks_attribs.push_back(t10);

      TaskAttributes t11;
      t11.advisory_total_num_threads = 8;
      t11.advisory_num_threads_per_group = 8;
      t11.task_type = 1;
      t11.buffer_binds.push_back({BufferEnum::Root,0});
      t11.buffer_binds.push_back({BufferEnum::GlobalTmps,1});
      t11.buffer_binds.push_back({BufferEnum::Context,2});
      t11.name = "substep_c4_0_k0006_vk_t01";
      t11.range_for_attribs = std::make_optional<TaskAttributes::RangeForAttributes>();
      t11.range_for_attribs->begin = 0;
      t11.range_for_attribs->end = 8;
      t11.range_for_attribs->const_begin = true;
      t11.range_for_attribs->const_end = true;
      reg_params1.kernel_attribs.tasks_attribs.push_back(t11);

      TaskAttributes t12;
      t12.advisory_total_num_threads = 256;
      t12.advisory_num_threads_per_group = 128;
      t12.task_type = 1;
      t12.buffer_binds.push_back({BufferEnum::Root,0});
      t12.buffer_binds.push_back({BufferEnum::GlobalTmps,1});
      t12.buffer_binds.push_back({BufferEnum::Context,2});
      t12.name = "substep_c4_0_k0006_vk_t02";
      t12.range_for_attribs = std::make_optional<TaskAttributes::RangeForAttributes>();
      t12.range_for_attribs->begin = 0;
      t12.range_for_attribs->end = 256;
      t12.range_for_attribs->const_begin = true;
      t12.range_for_attribs->const_end = true;
      reg_params1.kernel_attribs.tasks_attribs.push_back(t12);

      TaskAttributes t13;
      t13.advisory_total_num_threads = 8;
      t13.advisory_num_threads_per_group = 8;
      t13.task_type = 1;
      t13.buffer_binds.push_back({BufferEnum::Root,0});
      t13.buffer_binds.push_back({BufferEnum::GlobalTmps,1});
      t13.buffer_binds.push_back({BufferEnum::Context,2});
      t13.name = "substep_c4_0_k0006_vk_t03";
      t13.range_for_attribs = std::make_optional<TaskAttributes::RangeForAttributes>();
      t13.range_for_attribs->begin = 0;
      t13.range_for_attribs->end = 8;
      t13.range_for_attribs->const_begin = true;
      t13.range_for_attribs->const_end = true;
      reg_params1.kernel_attribs.tasks_attribs.push_back(t13);


      KernelContextAttributes::ArgAttributes arg1;
      arg1.stride = size_t(4);
      arg1.offset_in_mem = size_t(0);
      arg1.index = 0;
      arg1.is_array = false;
      reg_params1.kernel_attribs.ctx_attribs.arg_attribs_vec_.push_back(arg1);
      reg_params1.kernel_attribs.ctx_attribs.args_bytes_=4;
      reg_params1.kernel_attribs.ctx_attribs.rets_bytes_=0;
      reg_params1.kernel_attribs.ctx_attribs.extra_args_bytes_=256;
      register_params_.push_back(std::move(reg_params1));
    }

    KernelHandle register_taichi_kernel(const RegisterParams &reg_params) {
      CompiledTaichiKernel::Params params;
      params.ti_kernel_attribs = &(reg_params.kernel_attribs);
//      params.snode_descriptors = snode_descriptors_;
      params.device = managed_device_->device();
      params.root_buffer = root_buffer_.get();
      params.global_tmps_buffer = global_tmps_buffer_.get();
      params.host_visible_mem_pool = host_visible_memory_pool_.get();

      for (const auto &attribs: reg_params.kernel_attribs.tasks_attribs) {
//        const auto &glsl_src = reg_params.task_glsl_source_codes[i];
        const auto &task_name = attribs.name;
//        auto spv_bin = spv_compiler_.compile(glsl_src, task_name).value();
        auto spv_bin = spv_compiler_.load(task_name).value();
        // If we can reach here, we have succeeded. Otherwise
        // std::optional::value() would have killed us.
//        TI_TRACE("Successfully compiled GLSL -> SPIR-V for task={}\n{}",
//                 task_name, glsl_src);
        params.spirv_bins.push_back(std::move(spv_bin));
      }
      KernelHandle res;
      res.id_ = ti_kernels_.size();
      auto ctk = std::make_unique<CompiledTaichiKernel>(params);
      ti_kernels_.push_back(std::move(ctk));
      return res;
    }

    void launch_kernel(KernelHandle handle, Context *host_ctx) {
      auto *ti_kernel = ti_kernels_[handle.id_].get();
      auto ctx_blitter = HostDeviceContextBlitter::maybe_make(
              &ti_kernel->ti_kernel_attribs().ctx_attribs, host_ctx,
              host_result_buffer_, ti_kernel->ctx_buffer());
      if (ctx_blitter) {
        assert(ti_kernel->ctx_buffer() != nullptr);
//        TI_ASSERT(ti_kernel->ctx_buffer() != nullptr);
        ctx_blitter->host_to_device();
      }

      stream_->launch(ti_kernel->command_buffer());
      num_pending_kernels_ += ti_kernel->num_vk_pipelines();
      if (ctx_blitter) {
        synchronize();
        ctx_blitter->device_to_host();
      }
      synchronize();
    }

    void free_command_buffer(KernelHandle handle){
      auto *ti_kernel = ti_kernels_[handle.id_].get();
      auto cb = ti_kernel->command_buffer();
      auto dev = managed_device_->get_device();
      auto cp = managed_device_->get_command_pool();
      vkFreeCommandBuffers(dev, cp, 1, &cb);
    }

private:
    void init_vk_buffers() {
//  #pragma message("Vulkan buffers size hardcoded")
      root_buffer_ = dev_local_memory_pool_->alloc_and_bind(16 * 1024 * 1024);
      global_tmps_buffer_ = dev_local_memory_pool_->alloc_and_bind(1024 * 1024);

      // Need to zero fill the buffers, otherwise there could be NaN.
      ClearBufferCommandBuilder cmd_builder{stream_->device()};
      auto clear_cmd = cmd_builder.build(
              /*buffers=*/{root_buffer_->buffer(), global_tmps_buffer_->buffer()});
      stream_->launch(clear_cmd);
      stream_->synchronize();
    }

    void synchronize() {
      if (num_pending_kernels_ == 0) {
        return;
      }

//      TI_AUTO_PROF;
//      StopWatch sw;
      stream_->synchronize();
//      TI_DEBUG("running {} kernels took {} us", num_pending_kernels_,
//               sw.GetMicros());
      num_pending_kernels_ = 0;
    }

    void init_memory_pool(const Params &params) {
      LinearVkMemoryPool::Params mp_params;
      mp_params.physical_device = managed_device_->physical_device();
      mp_params.device = managed_device_->device()->device();
//    #pragma message("Vulkan memory pool size hardcoded to 64MB")
      mp_params.pool_size = 64 * 1024 * 1024;
      mp_params.required_properties = VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT;
      mp_params.compute_queue_family_index =
              managed_device_->queue_family_indices().compute_family.value();

      auto &buf_creation_template = mp_params.buffer_creation_template;
      buf_creation_template.sType = VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO;
      buf_creation_template.pNext = nullptr;
      buf_creation_template.flags = 0;
      buf_creation_template.size = 0;
      buf_creation_template.usage =
              (VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_SRC_BIT |
               VK_BUFFER_USAGE_TRANSFER_DST_BIT);
      buf_creation_template.sharingMode = VK_SHARING_MODE_EXCLUSIVE;
      buf_creation_template.queueFamilyIndexCount = 1;
      buf_creation_template.pQueueFamilyIndices = nullptr;
      dev_local_memory_pool_ = LinearVkMemoryPool::try_make(mp_params);
      assert(dev_local_memory_pool_);

      mp_params.required_properties = (VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT |
                                       VK_MEMORY_PROPERTY_HOST_COHERENT_BIT );
//      |
//                                       VK_MEMORY_PROPERTY_HOST_CACHED_BIT);

      buf_creation_template.usage =
              (VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_SRC_BIT);
      host_visible_memory_pool_ = LinearVkMemoryPool::try_make(mp_params);
      assert(host_visible_memory_pool_);
    }

    uint64_t *const host_result_buffer_;
    std::vector<std::unique_ptr<CompiledTaichiKernel>> ti_kernels_;
    std::vector<RegisterParams> register_params_;
    std::unique_ptr<ManagedVulkanDevice> managed_device_{nullptr};
    std::unique_ptr<VulkanStream> stream_{nullptr};
    GlslToSpirvCompiler spv_compiler_;

    std::unique_ptr<LinearVkMemoryPool> dev_local_memory_pool_;
    std::unique_ptr<VkBufferWithMemory> root_buffer_;
    std::unique_ptr<VkBufferWithMemory> global_tmps_buffer_;
    std::unique_ptr<LinearVkMemoryPool> host_visible_memory_pool_;
    int num_pending_kernels_{0};
};

using FunctionType = std::function<void(Context &)>;
using Result = RegisterParams;
//FunctionType compile_to_executable(Kernel *kernel, Test *runtime) {
//  const auto id = get_kernel_id();
//  const auto taichi_kernel_name(fmt::format("{}_k{:04d}_vk", kernel->name, id));
////  TI_INFO("VK codegen for Taichi kernel={}", taichi_kernel_name);
////  KernelCodegen::Params params;
////  params.ti_kernel_name = taichi_kernel_name;
////  params.kernel = kernel;
////  params.compiled_structs = compiled_structs;
////  KernelCodegen codegen(params);
////  auto res = codegen.run();
////  auto handle = runtime->register_taichi_kernel(std::move(res));
////  return [runtime, handle, taichi_kernel_name](Context &ctx) {
////    runtime->launch_kernel(handle, &ctx);
////  };
//}

bool InitVulkan(android_app* app) {
//  if (!gRenderer.Init(app, kTAG)) {
//    return false;
//  }
//
//  // create vertex / index buffer
//  const std::vector<float> vertexData = {
//    -1.0f, 1.0f, 0.0f,
//    1.0f, 1.0f, 0.0f,
//    0.0f, -1.0f, 0.0f,
//  };
//
  // 1. instantiate a surface
//  auto surf = std::make_shared<RenderSurface>();
//  surf->mVertexCount = 3;
//  surf->mInstanceCount = 1;
//  surf->mItemSize = 3;
//
//  // 2. create vertex / index buffer
//  gRenderer.CreateVertexBuffer(vertexData, surf);
//
//  // 3. create shader
//  gRenderer.CreateGraphicsPipeline("shaders/tri.vert.spv",
//          "shaders/tri.frag.spv", surf);
//
//  // 4. add surfs to the render list.
//  gRenderer.AddSurface(surf);
//
//  // 5. create command buffer and render passes.
//  gRenderer.ConstructRenderPass();

  Params params_;
  Test vk_app(params_, app);
  vk_app.run();
//  try {
////    vk_app.run();
//  } catch (const std::exception& e) {
//    std::cerr << e.what() << std::endl;
//    return EXIT_FAILURE;
//  }

  return EXIT_SUCCESS;
//  return true;
}

bool IsVulkanReady() {
  return true;
//  return gRenderer.IsReady();
}

void TerminateVulkan() {
//  gRenderer.Terminate();
}

bool VulkanRenderFrame() {
//  gRenderer.RenderFrame();
  return true;
}
#pragma once

#include <vector>
#include <unordered_map>
//#include <vulkan/vulkan.h>
//#include <vulkan/vulkan_core.h>
#include <vulkan_wrapper.h>

#include "vulkan_api.h"
#include "vulkan_common.h"
#include "vulkan_utils.h"
#include "kernel_utils.h"
#include "vulkan_simple_memory_pool.h"

using BufferEnum = TaskAttributes::Buffers;
using InputBuffersMap = std::unordered_map<BufferEnum, VkBufferWithMemory *>;

struct RegisterParams {
    TaichiKernelAttributes kernel_attribs;
    std::vector<std::string> task_spirv_source_codes;
    std::vector<std::string> task_glsl_source_codes;
};

// Info for launching a compiled Taichi kernel, which consists of a series of
// Vulkan pipelines.

class CompiledTaichiKernel {
public:
    struct Params {
        const TaichiKernelAttributes *ti_kernel_attribs{nullptr};
        std::vector<GlslToSpirvCompiler::SpirvBinary> spirv_bins;
//    const SNodeDescriptorsMap *snode_descriptors{nullptr};

        const VulkanDevice *device{nullptr};
        VkBufferWithMemory *root_buffer{nullptr};
        VkBufferWithMemory *global_tmps_buffer{nullptr};
        LinearVkMemoryPool *host_visible_mem_pool{nullptr};
    };

    CompiledTaichiKernel(const Params &ti_params):
            ti_kernel_attribs_(*ti_params.ti_kernel_attribs){
        InputBuffersMap input_buffers = {
                {BufferEnum::Root, ti_params.root_buffer},
                {BufferEnum::GlobalTmps, ti_params.global_tmps_buffer},
        };
        if (!ti_kernel_attribs_.ctx_attribs.empty()) {
            const auto ctx_sz = ti_kernel_attribs_.ctx_attribs.total_bytes();
            ctx_buffer_ = ti_params.host_visible_mem_pool->alloc_and_bind(ctx_sz);
            input_buffers[BufferEnum::Context] = ctx_buffer_.get();
        }

        const auto &task_attribs = ti_kernel_attribs_.tasks_attribs;
        const auto &spirv_bins = ti_params.spirv_bins;

        VulkanComputeCommandBuilder cmd_builder(ti_params.device);
        for (int i = 0; i < task_attribs.size(); ++i) {
            const auto &attribs = task_attribs[i];
            VulkanPipeline::Params vp_params;
            vp_params.device = ti_params.device;
            for (const auto &bb : task_attribs[i].buffer_binds) {
                auto b = (uint32_t)bb.binding;
                auto a = input_buffers.at(bb.type);
                auto c = a->buffer();

                VulkanPipeline::BufferBinding bb_{c, b};
                vp_params.buffer_bindings.push_back(bb_);
            }
            vp_params.code = SpirvCodeView(spirv_bins[i]);
            auto vp = std::make_unique<VulkanPipeline>(vp_params);
            const int group_x = attribs.advisory_total_num_threads /
                                attribs.advisory_num_threads_per_group;
            cmd_builder.append(*vp, group_x);
            vk_pipelines_.push_back(std::move(vp));
        }
        command_buffer_ = cmd_builder.build();
    }

    const TaichiKernelAttributes &ti_kernel_attribs() const {
        return ti_kernel_attribs_;
    }

    size_t num_vk_pipelines() const {
        return vk_pipelines_.size();
    }

    VkBufferWithMemory *ctx_buffer() const {
        return ctx_buffer_.get();
    }

    VkCommandBuffer command_buffer() const {
        return command_buffer_;
    }

private:
    TaichiKernelAttributes ti_kernel_attribs_;

    // Right now |ctx_buffer_| is allocated from a HOST_VISIBLE|COHERENT
    // memory, because we do not do computation on this buffer anyway, and it may
    // not worth the effort doing another hop via a staging buffer.
    // TODO: Provide an option to use staging buffer. This could be useful if the
    // kernel does lots of IO on the context buffer, e.g., copy a large np array.
    std::unique_ptr<VkBufferWithMemory> ctx_buffer_{nullptr};
    std::vector<std::unique_ptr<VulkanPipeline>> vk_pipelines_;

    // VkCommandBuffers are destroyed when the underlying command pool is
    // destroyed.
    // https://vulkan-tutorial.com/Drawing_a_triangle/Drawing/Command_buffers#page_Command-buffer-allocation
    VkCommandBuffer command_buffer_{VK_NULL_HANDLE};
};

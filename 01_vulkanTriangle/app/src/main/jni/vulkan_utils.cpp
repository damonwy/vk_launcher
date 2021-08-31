#include "vulkan_utils.h"
#include <iostream>
#include <fstream>
#include <iterator>
#include <stdio.h>
#include <cerrno>
#include <cstring>
#include <Logger.h>

//#include <spirv-tools/libspirv.hpp>

std::vector<VkExtensionProperties> GetInstanceExtensionProperties() {
    constexpr char *kNoLayerName = nullptr;
    uint32_t count = 0;
    vkEnumerateInstanceExtensionProperties(kNoLayerName, &count, nullptr);
    std::vector<VkExtensionProperties> extensions(count);
    vkEnumerateInstanceExtensionProperties(kNoLayerName, &count,
                                           extensions.data());
    return extensions;
}

std::vector<VkExtensionProperties> GetDeviceExtensionProperties(
        VkPhysicalDevice physicalDevice) {
    constexpr char *kNoLayerName = nullptr;
    uint32_t count = 0;
    vkEnumerateDeviceExtensionProperties(physicalDevice, kNoLayerName, &count,
                                         nullptr);
    std::vector<VkExtensionProperties> extensions(count);
    vkEnumerateDeviceExtensionProperties(physicalDevice, kNoLayerName, &count,
                                         extensions.data());
    return extensions;
}

GlslToSpirvCompiler::GlslToSpirvCompiler(){

}

std::optional<GlslToSpirvCompiler::SpirvBinary> GlslToSpirvCompiler::load(const std::string &task_name){
//    auto res = std::make_optional<GlslToSpirvCompiler::SpirvBinary>();
//    std::ifstream stream(task_name+".spv", std::ios::in|std::ios::binary);
//    std::vector<char> contents((std::istreambuf_iterator<char>(stream)), std::istreambuf_iterator<char>());

//    FILE* file = std::fopen(("/storage/emulated/0/Android/data/com.example.vulkanTriangle/files/shaders/"+task_name+".spv").c_str(), "rb");
    std::ifstream file(("/storage/emulated/0/Android/data/com.example.vulkanTriangle/files/shaders/"+task_name+".spv").c_str(), std::ios::binary);
    if (!file){
        LOG_E("CANNOT: {} {}", std::strerror(errno), "a");
    }

//    std::vector<char> spirv;
//    char buffer[1024];
//    int len;
//    while(len = std::fread(buffer, sizeof(char), sizeof(buffer), file) != EOF){
//        spirv.insert(spirv.end(), buffer, buffer+len);
//
//    }
//    std::fclose(file);

    if (!file.is_open())
        std::cout<<"ERROR in GlslToSpirvCompiler::load(): Could not load '" + std::string(task_name) + "' file." << std::endl;

    std::vector<char> spirv((std::istreambuf_iterator<char>(file)), (std::istreambuf_iterator<char>()));
    std::vector<uint32_t> spv(spirv.size() / sizeof(uint32_t));
//    // Copy data from the char-vector to a new uint32_t-vector
    memcpy(spv.data(), spirv.data(), spirv.size());

    return spv;
//
//    for(auto i: contents) {
//        int value = i;
//        std::cout << "data: " << value << std::endl;
//        if(j%4==0){
//            c0 = value;
//        }else if(j%4==1){
//            c1 = value;
//        }else if(j%4==2){
//            c2 = value;
//        }else if(j%4==3){
//            c3 = value;
//            uint32_t combined = (c3 << 24 | c2 << 16 | c1 << 8 | c0);
//            ret.push_back(combined);
//            std::cout << "combined: " << combined << std::endl;
//            std::cout << j << std::endl;
//        }
//        std::cout << j << std::endl;
//        j++;
//
//    }
//
//    std::cout << "file size: " << ret.size() << std::endl;
//    return ret;
}

//GlslToSpirvCompiler::GlslToSpirvCompiler(const ErrorHandler &err_handler)
//        : err_handler_(err_handler) {
//    opts_.SetTargetEnvironment(shaderc_target_env_vulkan,VulkanEnvSettings::kShadercEnvVersion());
//    opts_.SetOptimizationLevel(shaderc_optimization_level_performance);
//}

//std::optional<GlslToSpirvCompiler::SpirvBinary> GlslToSpirvCompiler::compile(
//        const std::string &glsl_src,
//        const std::string &shader_name) {
//    auto spv_result =
//            compiler_.CompileGlslToSpv(glsl_src, shaderc_glsl_default_compute_shader,
//                    /*input_file_name=*/shader_name.c_str(),
//                    /*entry_point_name=*/"main", opts_);
//    if (spv_result.GetCompilationStatus() != shaderc_compilation_status_success) {
//        err_handler_(glsl_src, shader_name, spv_result.GetErrorMessage());
//        return std::nullopt;
//    }
//    return SpirvBinary(spv_result.begin(), spv_result.end());
//}

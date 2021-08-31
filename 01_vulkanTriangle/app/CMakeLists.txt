cmake_minimum_required(VERSION 3.4.1)

# build native_app_glue as a static lib
set(APP_GLUE_DIR ${ANDROID_NDK}/sources/android/native_app_glue)
include_directories(${APP_GLUE_DIR})
add_library(app-glue STATIC ${APP_GLUE_DIR}/android_native_app_glue.c)

# build vulkan app
set(SRC_JNI_DIR src/main/jni)
set(SRC_RENDERER_DIR ../../common/renderer)
set(WRAPPER_DIR ../../common/vulkan_wrapper)
set(UTILS_DIR ../../common/utils)
set(THIRD_PARTY_DIR ../../third_party)


FILE(GLOB_RECURSE SOURCES ${SRC_JNI_DIR}/*.cpp)
FILE(GLOB_RECURSE HEADERS ${SRC_JNI_DIR}/*.h)

add_library(vkexamples SHARED
        ${SRC_JNI_DIR}/vulkan_api.cpp
        ${SRC_JNI_DIR}/kernel_utils.cpp
        ${SRC_JNI_DIR}/memory_pool.cpp
        ${SRC_JNI_DIR}/type.cpp
        ${SRC_JNI_DIR}/type_factory.cpp
        ${SRC_JNI_DIR}/type_utils.cpp
        ${SRC_JNI_DIR}/unified_allocator.cpp
        ${SRC_JNI_DIR}/vulkan_simple_memory_pool.cpp
        ${SRC_JNI_DIR}/vulkan_utils.cpp

        ${SRC_JNI_DIR}/AndroidMain.cpp
            ${SRC_JNI_DIR}/VulkanMain.cpp
            ${SRC_RENDERER_DIR}/VulkanRenderer.cpp
            ${WRAPPER_DIR}/vulkan_wrapper.cpp
            ${UTILS_DIR}/Platform.cpp)

include_directories(${WRAPPER_DIR}
                    ${UTILS_DIR}
                    ${SRC_RENDERER_DIR}
        ${SRC_JNI_DIR}/vulkan_api.h
        ${SRC_JNI_DIR}/kernel_utils.h
        ${SRC_JNI_DIR}/memory_pool.h
        ${SRC_JNI_DIR}/type.h
        ${SRC_JNI_DIR}/type_factory.h
        ${SRC_JNI_DIR}/type_utils.h
        ${SRC_JNI_DIR}/unified_allocator.h
        ${SRC_JNI_DIR}/vulkan_simple_memory_pool.h
        ${SRC_JNI_DIR}/vulkan_utils.h
        ${SRC_JNI_DIR}/arithmetic.h
        ${SRC_JNI_DIR}/common.h
        ${SRC_JNI_DIR}/constants.h
        ${SRC_JNI_DIR}/context.h
        ${SRC_JNI_DIR}/data_type_inc.h
        ${SRC_JNI_DIR}/data_type_utils.h
        ${SRC_JNI_DIR}/mem_request.h
        ${SRC_JNI_DIR}/runtime.h
        ${SRC_JNI_DIR}/virtual_memory.h


        ${THIRD_PARTY_DIR}/gfx-math/include
                    ${THIRD_PARTY_DIR}/KTX-Software/include
                    ${THIRD_PARTY_DIR}/spdlog/include)

# Add third party libraries
add_subdirectory(../../third_party third_party)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -Wall \
                     -DVK_USE_PLATFORM_ANDROID_KHR")
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -u ANativeActivity_onCreate")

# 'ktx' is built from third_party.
target_link_libraries(vkexamples app-glue log android ktx)

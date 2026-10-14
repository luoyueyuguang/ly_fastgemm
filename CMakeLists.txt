cmake_minimum_required(VERSION 3.14)
project(lygemm VERSION 2.0.0 LANGUAGES CXX)

# Header only: the interface target is the whole package.
add_library(lygemm INTERFACE)
add_library(lygemm::lygemm ALIAS lygemm)
target_include_directories(lygemm INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)
target_compile_features(lygemm INTERFACE cxx_std_17)

# The instruction set is the consumer's choice: no vector flags means the
# portable scalar backend, -mavx2 -mfma the 256 bit one, -mavx512f -mavx512dq
# the 512 bit one.  Note that AVX2 usually wins on shapes narrower than eight
# columns even on AVX-512 hardware, see the README.

include(GNUInstallDirs)
install(TARGETS lygemm EXPORT lygemmTargets)
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT lygemmTargets NAMESPACE lygemm::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/lygemm)

option(LYGEMM_BUILD_TESTS "Build the lygemm correctness suite" ON)
if(LYGEMM_BUILD_TESTS)
  enable_testing()
  add_executable(lygemm_test_correctness tests/test_correctness.cpp)
  target_link_libraries(lygemm_test_correctness PRIVATE lygemm)
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(lygemm_test_correctness PRIVATE
      -O3 -Wall -Wextra $<$<BOOL:${LYGEMM_TEST_ISA}>:${LYGEMM_TEST_ISA}>)
    if(NOT LYGEMM_TEST_ISA)
      target_compile_options(lygemm_test_correctness PRIVATE -mavx2 -mfma)
    endif()
  endif()
  add_test(NAME correctness COMMAND lygemm_test_correctness)
endif()

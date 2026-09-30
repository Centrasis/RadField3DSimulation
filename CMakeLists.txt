cmake_minimum_required(VERSION 3.24...3.27)
project(RadField3DSimulation LANGUAGES CXX)


SET(CMAKE_CONFIGURATION_TYPES "Release;RelWithDebInfo")
SET(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(FetchContent)
set(FETCHCONTENT_QUIET TRUE)
# Some fetched dependencies (e.g. nlohmann/json 3.11.2) declare cmake_minimum_required below 3.5, which CMake >= 4
# rejects unless a policy floor is given.
if(NOT DEFINED CMAKE_POLICY_VERSION_MINIMUM)
  set(CMAKE_POLICY_VERSION_MINIMUM 3.5)
endif()

option(BUILD_TESTS "Enable Tests using Google Test" OFF)
set(PROJECT_VERSION 1.0.0)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# ---------------------------------------------------------------------------
# Sanitizer build (for hunting the non-deterministic MT segfault).
#   cmake -S . -B build-asan -DCMAKE_BUILD_TYPE=RelWithDebInfo -DRS_SANITIZER=address
#   cmake -S . -B build-tsan -DCMAKE_BUILD_TYPE=RelWithDebInfo -DRS_SANITIZER=thread
# address = use-after-free / heap overflow / double-free (start here — it targets the
#           memory corruption that produces the SIGSEGV); thread = data races.
# Set here, before any target or FetchContent subdirectory is defined, so every TU
# (ours + radfiled3d, where the threaded scoring lives) is instrumented and the
# runtime is linked into the RadField3D executable. Geant4 itself is a prebuilt lib and
# stays uninstrumented — races inside our code are still caught; ASan/TSan tolerate the gap.
# Runtime hints: ASan  -> ASAN_OPTIONS=detect_leaks=0:abort_on_error=1  (Geant4 keeps
#                         singletons alive at exit, which LeakSanitizer would flag)
#                TSan  -> TSAN_OPTIONS=suppressions=<file>:history_size=7  to mute Geant4-internal races
# ---------------------------------------------------------------------------
set(RS_SANITIZER "OFF" CACHE STRING "Build with a sanitizer: OFF, address, thread, undefined, or address,undefined")
set_property(CACHE RS_SANITIZER PROPERTY STRINGS OFF address thread undefined address,undefined)
if(NOT RS_SANITIZER STREQUAL "OFF")
  if(MSVC)
    message(FATAL_ERROR "RS_SANITIZER is wired up for GCC/Clang only (the Linux cluster build).")
  endif()
  if(RS_SANITIZER MATCHES "address" AND RS_SANITIZER MATCHES "thread")
    message(FATAL_ERROR "AddressSanitizer and ThreadSanitizer are mutually exclusive; build them in separate trees.")
  endif()
  message(STATUS "Building with sanitizer: ${RS_SANITIZER}")
  # -UNDEBUG re-enables the detector's assert()s under RelWithDebInfo/Release (which define NDEBUG);
  # -fno-omit-frame-pointer / -fno-optimize-sibling-calls keep the crash stacks readable.
  add_compile_options(-fsanitize=${RS_SANITIZER} -fno-omit-frame-pointer -fno-optimize-sibling-calls -g -UNDEBUG)
  add_link_options(-fsanitize=${RS_SANITIZER})
endif()

message(STATUS "Download Eigen...")
FetchContent_Declare(
  eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen
  GIT_TAG        3.4.0
)
FetchContent_GetProperties(eigen)
if(NOT eigen_POPULATED)
  FetchContent_Populate(eigen)
  add_subdirectory(${eigen_SOURCE_DIR} ${eigen_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()

message(STATUS "Download Assimp...")
FetchContent_Declare(
  assimp
  GIT_REPOSITORY https://github.com/assimp/assimp.git
  GIT_TAG        v5.2.5
)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
set(ASSIMP_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(ASSIMP_INJECT_DEBUG_POSTFIX OFF CACHE BOOL "" FORCE)
set(ASSIMP_INSTALL OFF CACHE BOOL "" FORCE)
# Only OBJ/STL phantom geometry is loaded. Building just those importers avoids compiling assimp's other
# format parsers (e.g. the MMD/pmx parser, which fails on modern GCC without <cstdint>).
set(ASSIMP_BUILD_ALL_IMPORTERS_BY_DEFAULT OFF CACHE BOOL "" FORCE)
set(ASSIMP_BUILD_ALL_EXPORTERS_BY_DEFAULT OFF CACHE BOOL "" FORCE)
set(ASSIMP_BUILD_OBJ_IMPORTER ON CACHE BOOL "" FORCE)
set(ASSIMP_BUILD_STL_IMPORTER ON CACHE BOOL "" FORCE)
FetchContent_GetProperties(assimp)
if(NOT assimp_POPULATED)
  FetchContent_Populate(assimp)
  add_subdirectory(${assimp_SOURCE_DIR} ${assimp_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()

if (MSVC)
  target_compile_options(assimp PRIVATE /W3)
else()
  target_compile_options(assimp PRIVATE -Wall -Wno-error)
endif()

message(STATUS "Download JSON...")
FetchContent_Declare(
  json
  GIT_REPOSITORY https://github.com/nlohmann/json.git
  GIT_TAG        v3.11.2
)
FetchContent_GetProperties(json)
if(NOT json_POPULATED)
  FetchContent_Populate(json)
  add_subdirectory(${json_SOURCE_DIR} ${json_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()

# --- OpenVDB: voxelizes the scene geometry into the stored field ------------------------------------
# An installed OpenVDB (FindOpenVDB.cmake on CMAKE_MODULE_PATH) and oneTBB (TBBConfig.cmake) are used when
# found; otherwise both are fetched and built statically (CPU-only core library, no optional codecs),
# so a fresh clone builds without any extra system packages.
find_package(TBB CONFIG QUIET COMPONENTS tbb)
if(NOT TBB_FOUND)
  message(STATUS "Download oneTBB...")
  set(TBB_TEST OFF CACHE BOOL "" FORCE)
  set(TBB_EXAMPLES OFF CACHE BOOL "" FORCE)
  set(TBB_STRICT OFF CACHE BOOL "" FORCE)
  set(TBB_INSTALL OFF CACHE BOOL "" FORCE)
  set(TBBMALLOC_BUILD OFF CACHE BOOL "" FORCE)
  FetchContent_Declare(
    TBB
    GIT_REPOSITORY https://github.com/uxlfoundation/oneTBB.git
    GIT_TAG        v2022.3.0
    GIT_SHALLOW    TRUE
    OVERRIDE_FIND_PACKAGE
  )
  FetchContent_MakeAvailable(TBB)
endif()

find_package(OpenVDB QUIET)
if(TARGET OpenVDB::openvdb)
  set(RS_OPENVDB_TARGET OpenVDB::openvdb)
else()
  message(STATUS "Download OpenVDB...")
  set(OPENVDB_BUILD_CORE ON CACHE BOOL "" FORCE)
  set(OPENVDB_CORE_SHARED OFF CACHE BOOL "" FORCE)
  set(OPENVDB_CORE_STATIC ON CACHE BOOL "" FORCE)
  set(OPENVDB_BUILD_BINARIES OFF CACHE BOOL "" FORCE)
  set(OPENVDB_BUILD_UNITTESTS OFF CACHE BOOL "" FORCE)
  set(OPENVDB_BUILD_DOCS OFF CACHE BOOL "" FORCE)
  set(OPENVDB_BUILD_PYTHON_MODULE OFF CACHE BOOL "" FORCE)
  set(OPENVDB_BUILD_NANOVDB OFF CACHE BOOL "" FORCE)
  set(OPENVDB_BUILD_AX OFF CACHE BOOL "" FORCE)
  set(OPENVDB_USE_DELAYED_LOADING OFF CACHE BOOL "" FORCE)
  set(OPENVDB_ENABLE_UNINSTALL OFF CACHE BOOL "" FORCE)
  set(OPENVDB_INSTALL_CMAKE_MODULES OFF CACHE BOOL "" FORCE)
  set(OPENVDB_FUTURE_DEPRECATION OFF CACHE BOOL "" FORCE)
  set(USE_BLOSC OFF CACHE BOOL "" FORCE)
  set(USE_ZLIB OFF CACHE BOOL "" FORCE)
  set(USE_CCACHE OFF CACHE BOOL "" FORCE)
  set(USE_EXPLICIT_INSTANTIATION OFF CACHE BOOL "" FORCE)
  set(CONCURRENT_MALLOC "None" CACHE STRING "" FORCE)
  FetchContent_Declare(
    openvdb
    GIT_REPOSITORY https://github.com/AcademySoftwareFoundation/openvdb.git
    GIT_TAG        v13.1.0
    GIT_SHALLOW    TRUE
    PATCH_COMMAND  ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/PatchOpenVDB.cmake
  )
  FetchContent_MakeAvailable(openvdb)
  set(RS_OPENVDB_TARGET openvdb_static)
endif()

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

option(WITH_GEANT4_UIVIS "Build with Geant4 UI and Vis drivers" OFF)

if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()


message(STATUS "Write dependencies file: ${PROJECT_SOURCE_DIR}/Dependencies.conf")
file(WRITE "${PROJECT_SOURCE_DIR}/Dependencies.conf" "${Geant4_INCLUDE_DIRS}/../../bin\n")
file(APPEND "${PROJECT_SOURCE_DIR}/Dependencies.conf" "bin\n")

set(RS_Enable_Linetracing, ON)
# Develop against a RadFiled3D checkout instead of the pinned submodule, e.g.
# -DRS_RADFILED3D_SOURCE_DIR=/path/to/RadFiled3D (in a separate build directory).
set(RS_RADFILED3D_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/targets/radiationfieldstorage" CACHE PATH "RadFiled3D sources to build against")
if(NOT EXISTS "${RS_RADFILED3D_SOURCE_DIR}/include/radfiled3d/voxel.hpp")
  message(FATAL_ERROR "RadFiled3D in ${RS_RADFILED3D_SOURCE_DIR} is older than 1.4.0 (no include/radfiled3d/). "
    "Update the submodule (git submodule update --remote targets/radiationfieldstorage) or point RS_RADFILED3D_SOURCE_DIR at a 1.4.0 checkout.")
endif()
add_subdirectory(${RS_RADFILED3D_SOURCE_DIR} ${CMAKE_BINARY_DIR}/targets/radiationfieldstorage)

include(${Geant4_USE_FILE})

set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

file(GLOB_RECURSE sources PROJECT_SOURCE_FILES ${PROJECT_SOURCE_DIR}/src/*.cpp ${PROJECT_SOURCE_DIR}/include/*.hpp)

set_target_properties(radfiled3d PROPERTIES FOLDER "Library")
set_target_properties(glm PROPERTIES FOLDER "ThirdParty")
if (TARGET zlibstatic)
  set_target_properties(zlibstatic PROPERTIES FOLDER "ThirdParty")
endif()

SET(src_list ${sources})
list(APPEND src_list ${include})
foreach(_source IN ITEMS ${src_list})
    get_filename_component(_source_path ${_source} PATH)
    file(RELATIVE_PATH _source_path_rel "${PROJECT_SOURCE_DIR}" "${_source_path}")
    string(REPLACE "/" "\\" _group_path "${_source_path_rel}")
    source_group(${_group_path} FILES ${_source})
endforeach()

add_library(libRadField3D STATIC ${sources})
target_include_directories(libRadField3D PUBLIC include)
if (MSVC)
  target_link_libraries(libRadField3D ${Geant4_LIBRARIES} glm assimp nlohmann_json::nlohmann_json Eigen3::Eigen radfiled3d ${RS_OPENVDB_TARGET})
else()
  target_link_libraries(libRadField3D ${Geant4_LIBRARIES} glm assimp nlohmann_json::nlohmann_json Eigen3::Eigen radfiled3d ${RS_OPENVDB_TARGET} stdc++fs)
endif()
if(WITH_GEANT4_UIVIS)
  target_compile_definitions(libRadField3D PUBLIC WITH_GEANT4_UIVIS)
endif()
set_target_properties(libRadField3D
  PROPERTIES
  ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/lib"
  LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/lib"
  RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin"
)
set_target_properties(libRadField3D PROPERTIES FOLDER "Library")

add_executable(RadField3D ${PROJECT_SOURCE_DIR}/targets/radfield3d/RadField3D.cpp ${sources})
target_include_directories(RadField3D PUBLIC include)
if (MSVC)
  target_link_libraries(RadField3D ${Geant4_LIBRARIES} glm assimp nlohmann_json::nlohmann_json Eigen3::Eigen radfiled3d ${RS_OPENVDB_TARGET} libRadField3D)
else()
  target_link_libraries(RadField3D ${Geant4_LIBRARIES} glm assimp nlohmann_json::nlohmann_json Eigen3::Eigen radfiled3d ${RS_OPENVDB_TARGET} libRadField3D stdc++fs)
endif()
if(WITH_GEANT4_UIVIS)
  target_compile_definitions(RadField3D PUBLIC WITH_GEANT4_UIVIS)
endif()
set_target_properties(RadField3D
  PROPERTIES
  ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/lib"
  LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/lib"
  RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin"
)

# --- RadField3D-Analytical: GPU analytical (non-Monte-Carlo) simulator ----------------------------
# A drop-in --binary for tools/create_dataset.py that generates low-accuracy fields fast, to
# pre-flight network trainability. Geant4-free (reuses only GeometryLoader + the mesh types); needs
# CUDA + cuFFT. Disabled by default; enable with -DBUILD_ANALYTICAL=ON.
option(BUILD_ANALYTICAL "Build the analytical CUDA ray-tracing simulator target" OFF)
if(BUILD_ANALYTICAL)
  enable_language(CUDA)
  if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES native)
  endif()
  find_package(CUDAToolkit REQUIRED)
  add_executable(RadField3D-Analytical
    ${PROJECT_SOURCE_DIR}/targets/analytical/main.cpp
    ${PROJECT_SOURCE_DIR}/targets/analytical/AnalyticalVoxelizer.cpp
    ${PROJECT_SOURCE_DIR}/targets/analytical/AnalyticalTracer.cu
    ${PROJECT_SOURCE_DIR}/src/GeometryLoader.cpp
    ${PROJECT_SOURCE_DIR}/src/Geometry.cpp
    ${PROJECT_SOURCE_DIR}/src/Collisions.cpp
    ${PROJECT_SOURCE_DIR}/src/RadiationSource.cpp
    ${PROJECT_SOURCE_DIR}/src/ProbabilityFunctions.cpp
    ${PROJECT_SOURCE_DIR}/src/Voxelization.cpp
  )
  target_include_directories(RadField3D-Analytical PRIVATE include)
  target_link_libraries(RadField3D-Analytical glm assimp nlohmann_json::nlohmann_json Eigen3::Eigen radfiled3d ${RS_OPENVDB_TARGET} CUDA::cufft)
  if(NOT MSVC)
    target_link_libraries(RadField3D-Analytical stdc++fs)
  endif()
  set_target_properties(RadField3D-Analytical PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin"
    CUDA_SEPARABLE_COMPILATION ON)
  set_target_properties(RadField3D-Analytical PROPERTIES FOLDER "Targets")
endif()

if(UNIX)
  set(CPACK_PACKAGE_NAME "RadField3D")
  set(CPACK_PACKAGE_VERSION ${PROJECT_VERSION})
  set(CPACK_PACKAGE_CONTACT "Felix Lehner <felix.lehner@ptb.de>")
  set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "RadField3D Simulation Package for simulating rotating x-ray tubes.")
  set(CPACK_PACKAGE_VENDOR "PTB")
  set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-${CMAKE_SYSTEM_PROCESSOR}")
  set(CPACK_VERBATIM_VARIABLES YES)
  set(CPACK_PACKAGE_INSTALL_DIRECTORY ${CPACK_PACKAGE_NAME})
  set(CPACK_PACKAGING_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
  SET(CPACK_OUTPUT_FILE_PREFIX "${CMAKE_SOURCE_DIR}/_packages")
  set(CPACK_STRIP_FILES YES)
  set(
    CPACK_INSTALL_DEFAULT_DIRECTORY_PERMISSIONS
      OWNER_READ OWNER_WRITE OWNER_EXECUTE
      GROUP_READ GROUP_EXECUTE
      WORLD_READ WORLD_EXECUTE
  )
  set(CPACK_COMPONENTS_GROUPING ALL_COMPONENTS_IN_ONE)
  set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS YES)
  set(CPACK_DEB_COMPONENT_INSTALL YES)

  if(EXISTS "/etc/debian_version")
    set(CPACK_DEBIAN_FILE_NAME DEB-DEFAULT)
    set(CPACK_GENERATOR "DEB")
    set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Felix Lehner")
    set(CPACK_DEBIAN_PACKAGE_DEPENDS "libc6 (>= 2.27), geant4, zlib1g, libexpat1, libexpat1-dev")
  elseif(EXISTS "/etc/redhat-release")
    set(CPACK_GENERATOR "RPM")
    set(CPACK_RPM_RELOCATABLE TRUE)
    set(CPACK_RPM_PACKAGE_RELOCATABLE_PATHS "/usr/bin/RadField3D")
    set(CPACK_RPM_PACKAGE_REQUIRES "libc >= 2.27, zlib, expat, expat-devel")
  endif()

  set(CPACK_PACKAGE_EXECUTABLES "RadField3D" "RadField3D")

  install(TARGETS RadField3D
    RUNTIME DESTINATION /usr/bin
  )

  include(CPack)
endif()

if(BUILD_TESTS)
  include(FetchContent)
  FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
  )

  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(googletest)

  enable_testing()

  function(link_test_target target_name)
    target_include_directories(${target_name} PUBLIC include)
    # The FetchContent gtest headers must win over any gtest in the Geant4 install prefix (e.g. a
    # conda env) — mixed gtest header/lib versions fail at link (MakeAndRegisterTestInfo).
    target_include_directories(${target_name} BEFORE PRIVATE ${googletest_SOURCE_DIR}/googletest/include)
    if (MSVC)
      target_link_libraries(${target_name} ${Geant4_LIBRARIES} glm assimp nlohmann_json::nlohmann_json Eigen3::Eigen radfiled3d ${RS_OPENVDB_TARGET} libRadField3D GTest::gtest_main)
    else()
      target_link_libraries(${target_name} ${Geant4_LIBRARIES} glm assimp nlohmann_json::nlohmann_json Eigen3::Eigen radfiled3d ${RS_OPENVDB_TARGET} libRadField3D GTest::gtest_main stdc++fs)
    endif()
    set_property(TARGET ${target_name} PROPERTY CXX_STANDARD ${CMAKE_CXX_STANDARD})
    set_target_properties(${target_name}
      PROPERTIES
      ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/lib"
      LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/lib"
      RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin"
    )
  endfunction()

  add_executable(WorldTest ${PROJECT_SOURCE_DIR}/targets/test/WorldTest.cpp)
  link_test_target(WorldTest)

  add_executable(CollisionTest ${PROJECT_SOURCE_DIR}/targets/test/CollisionTest.cpp)
  link_test_target(CollisionTest)

  add_executable(HistogramStatisticsTest ${PROJECT_SOURCE_DIR}/targets/test/HistogramStatisticsTest.cpp)
  link_test_target(HistogramStatisticsTest)

  add_executable(SourceDistributionTests ${PROJECT_SOURCE_DIR}/targets/test/SourceDistributionTests.cpp)
  link_test_target(SourceDistributionTests)

  add_executable(DetectorTest ${PROJECT_SOURCE_DIR}/targets/test/DetectorTest.cpp)
  link_test_target(DetectorTest)

  add_executable(VoxelizationTest ${PROJECT_SOURCE_DIR}/targets/test/VoxelizationTest.cpp)
  link_test_target(VoxelizationTest)

  add_executable(GeometryChannelTest ${PROJECT_SOURCE_DIR}/targets/test/GeometryChannelTest.cpp)
  link_test_target(GeometryChannelTest)

  add_executable(PatientScoringTest ${PROJECT_SOURCE_DIR}/targets/test/PatientScoringTest.cpp)
  link_test_target(PatientScoringTest)

  add_executable(GeometryLoaderTest ${PROJECT_SOURCE_DIR}/targets/test/GeometryLoaderTest.cpp)
  link_test_target(GeometryLoaderTest)

  add_executable(FieldAppendTest ${PROJECT_SOURCE_DIR}/targets/test/FieldAppendTest.cpp)
  link_test_target(FieldAppendTest)

  add_executable(SpectrumSamplingTest ${PROJECT_SOURCE_DIR}/targets/test/SpectrumSamplingTest.cpp)
  link_test_target(SpectrumSamplingTest)

  add_executable(PhysicsListTest ${PROJECT_SOURCE_DIR}/targets/test/PhysicsListTest.cpp)
  link_test_target(PhysicsListTest)

  add_executable(VMFTrainerTest ${PROJECT_SOURCE_DIR}/targets/test/VMFTrainerTest.cpp)
  link_test_target(VMFTrainerTest)

  set_target_properties(WorldTest PROPERTIES FOLDER "Tests")
  set_target_properties(CollisionTest PROPERTIES FOLDER "Tests")
  set_target_properties(DetectorTest PROPERTIES FOLDER "Tests")
  set_target_properties(VoxelizationTest PROPERTIES FOLDER "Tests")
  set_target_properties(GeometryChannelTest PROPERTIES FOLDER "Tests")
  set_target_properties(PatientScoringTest PROPERTIES FOLDER "Tests")
  set_target_properties(GeometryLoaderTest PROPERTIES FOLDER "Tests")
  set_target_properties(FieldAppendTest PROPERTIES FOLDER "Tests")
  set_target_properties(SpectrumSamplingTest PROPERTIES FOLDER "Tests")
  set_target_properties(PhysicsListTest PROPERTIES FOLDER "Tests")
  set_target_properties(VMFTrainerTest PROPERTIES FOLDER "Tests")
  set_target_properties(HistogramStatisticsTest PROPERTIES FOLDER "Tests")
  set_target_properties(SourceDistributionTests PROPERTIES FOLDER "Tests")


  set_target_properties(gtest PROPERTIES FOLDER "ThirdParty/GTest")
  set_target_properties(gtest_main PROPERTIES FOLDER "ThirdParty/GTest")
  set_target_properties(gmock PROPERTIES FOLDER "ThirdParty/GTest")
  set_target_properties(gmock_main PROPERTIES FOLDER "ThirdParty/GTest")


  include(GoogleTest)

  gtest_discover_tests(CollisionTest)
  gtest_discover_tests(VoxelizationTest)
  gtest_discover_tests(GeometryLoaderTest)
  gtest_discover_tests(FieldAppendTest)
  gtest_discover_tests(SpectrumSamplingTest)
  gtest_discover_tests(VMFTrainerTest)
 # gtest_discover_tests(RSLibLoadTest)
endif()

set_target_properties(assimp PROPERTIES FOLDER "ThirdParty")
cmake_minimum_required(VERSION 3.16...3.27)
project(RadField3DSimulation LANGUAGES CXX)


SET(CMAKE_CONFIGURATION_TYPES "Release;RelWithDebInfo;Debug")
SET(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(FetchContent)
set(FETCHCONTENT_QUIET TRUE)

option(BUILD_TESTS "Enable Tests using Google Test" OFF)
set(PROJECT_VERSION 1.0.0)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

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
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/targets/radiationfieldstorage)

include(${Geant4_USE_FILE})

set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

file(GLOB_RECURSE sources PROJECT_SOURCE_FILES ${PROJECT_SOURCE_DIR}/src/*.cpp ${PROJECT_SOURCE_DIR}/include/*.hpp)

set_target_properties(libRadFiled3D PROPERTIES FOLDER "Library")
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

add_executable(RadField3D ${sources})
target_include_directories(RadField3D PUBLIC include PRIVATE src)
if (MSVC)
  target_link_libraries(RadField3D ${Geant4_LIBRARIES} glm assimp nlohmann_json::nlohmann_json Eigen3::Eigen libRadFiled3D)
else()
  target_link_libraries(RadField3D ${Geant4_LIBRARIES} glm assimp nlohmann_json::nlohmann_json Eigen3::Eigen libRadFiled3D stdc++fs)
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

  add_executable(WorldTest ${PROJECT_SOURCE_DIR}/targets/test/WorldTest.cpp)
  target_include_directories(WorldTest PUBLIC include PRIVATE src)
  if (MSVC)
    target_link_libraries(WorldTest ${Geant4_LIBRARIES} glm assimp nlohmann_json::nlohmann_json Eigen3::Eigen libRadFiled3D)
  else()
    target_link_libraries(WorldTest ${Geant4_LIBRARIES} glm assimp nlohmann_json::nlohmann_json Eigen3::Eigen libRadFiled3D stdc++fs)
  endif()
  set_property(TARGET WorldTest PROPERTY CXX_STANDARD ${CMAKE_CXX_STANDARD})
  set_target_properties(WorldTest
      PROPERTIES
      ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/lib"
      LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/lib"
      RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin"
  )

  add_executable(CollisionTest ${PROJECT_SOURCE_DIR}/targets/test/CollisionTest.cpp)
  target_include_directories(CollisionTest PUBLIC include PRIVATE src)
  if (MSVC)
    target_link_libraries(CollisionTest ${Geant4_LIBRARIES} glm assimp nlohmann_json::nlohmann_json Eigen3::Eigen libRadFiled3D GTest::gtest_main)
  else()
    target_link_libraries(CollisionTest ${Geant4_LIBRARIES} glm assimp nlohmann_json::nlohmann_json Eigen3::Eigen libRadFiled3D GTest::gtest_main stdc++fs)
  endif()
  set_property(TARGET CollisionTest PROPERTY CXX_STANDARD ${CMAKE_CXX_STANDARD})
  set_target_properties(CollisionTest
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/lib"
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin"
  )

  add_executable(HistogramStatisticsTest ${PROJECT_SOURCE_DIR}/targets/test/HistogramStatisticsTest.cpp)
  target_include_directories(HistogramStatisticsTest PUBLIC include PRIVATE src)
  if (MSVC)
    target_link_libraries(HistogramStatisticsTest ${Geant4_LIBRARIES} glm assimp nlohmann_json::nlohmann_json Eigen3::Eigen libRadFiled3D GTest::gtest_main)
  else()
    target_link_libraries(HistogramStatisticsTest ${Geant4_LIBRARIES} glm assimp nlohmann_json::nlohmann_json Eigen3::Eigen libRadFiled3D GTest::gtest_main stdc++fs)
  endif()
  set_property(TARGET HistogramStatisticsTest PROPERTY CXX_STANDARD ${CMAKE_CXX_STANDARD})
  set_target_properties(HistogramStatisticsTest
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/lib"
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin"
  )

  add_executable(RSLibPCADirectionTest ${PROJECT_SOURCE_DIR}/targets/test/PCADirectionTests.cpp)
  target_include_directories(RSLibPCADirectionTest PUBLIC include PRIVATE src)
  if (MSVC)
    target_link_libraries(RSLibPCADirectionTest ${Geant4_LIBRARIES} glm assimp nlohmann_json::nlohmann_json Eigen3::Eigen libRadFiled3D GTest::gtest_main)
  else()
    target_link_libraries(RSLibPCADirectionTest ${Geant4_LIBRARIES} glm assimp nlohmann_json::nlohmann_json Eigen3::Eigen libRadFiled3D GTest::gtest_main stdc++fs)
  endif()
  set_property(TARGET RSLibPCADirectionTest PROPERTY CXX_STANDARD ${CMAKE_CXX_STANDARD})
  set_target_properties(RSLibPCADirectionTest
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/lib"
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin"
  )

  add_executable(RSLibLoadTest ${PROJECT_SOURCE_DIR}/targets/test/LoadTest.cpp)
  target_include_directories(RSLibLoadTest PUBLIC include PRIVATE src)
  if (MSVC)
    target_link_libraries(RSLibLoadTest ${Geant4_LIBRARIES} glm assimp nlohmann_json::nlohmann_json Eigen3::Eigen libRadFiled3D GTest::gtest_main)
  else()
    target_link_libraries(RSLibLoadTest ${Geant4_LIBRARIES} glm assimp nlohmann_json::nlohmann_json Eigen3::Eigen libRadFiled3D GTest::gtest_main stdc++fs)
  endif()
  set_property(TARGET RSLibLoadTest PROPERTY CXX_STANDARD ${CMAKE_CXX_STANDARD})
  set_target_properties(RSLibLoadTest
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/lib"
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin"
  )

  set_target_properties(WorldTest PROPERTIES FOLDER "Tests")
  set_target_properties(CollisionTest PROPERTIES FOLDER "Tests")
  set_target_properties(RSLibLoadTest PROPERTIES FOLDER "Tests")
  set_target_properties(RSLibPCADirectionTest PROPERTIES FOLDER "Tests")
  set_target_properties(HistogramStatisticsTest PROPERTIES FOLDER "Tests")


  set_target_properties(gtest PROPERTIES FOLDER "ThirdParty/GTest")
  set_target_properties(gtest_main PROPERTIES FOLDER "ThirdParty/GTest")
  set_target_properties(gmock PROPERTIES FOLDER "ThirdParty/GTest")
  set_target_properties(gmock_main PROPERTIES FOLDER "ThirdParty/GTest")


  include(GoogleTest)

  gtest_discover_tests(CollisionTest)
 # gtest_discover_tests(RSLibLoadTest)
endif()

set_target_properties(assimp PROPERTIES FOLDER "ThirdParty")
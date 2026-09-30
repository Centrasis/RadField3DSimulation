# Lets OpenVDB build as a FetchContent subproject: its core CMakeLists resolves its own files through
# ${CMAKE_SOURCE_DIR}, which only points at OpenVDB when it is the top-level project. Runs as the
# FetchContent PATCH_COMMAND with the OpenVDB source tree as working directory; safe to re-run.
foreach(file IN ITEMS openvdb/openvdb/CMakeLists.txt cmake/OpenVDBUtils.cmake)
  file(READ "${file}" content)
  string(REPLACE "\${CMAKE_SOURCE_DIR}" "\${OpenVDB_SOURCE_DIR}" patched "${content}")
  if(NOT patched STREQUAL content)
    file(WRITE "${file}" "${patched}")
  endif()
endforeach()

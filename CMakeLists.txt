cmake_minimum_required(VERSION 3.20)
project(oILAB VERSION 1.0 LANGUAGES C CXX)

# ---------- C++ standard ----------
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release)
endif()

if(MSVC)
    add_compile_options(/bigobj)
endif()

# ---------- Dependencies ----------
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# ---------- Library type ----------
option(BUILD_SHARED_LIBS "Build oILAB as a shared library" ON)

# ---------- Subdirectories ----------
add_subdirectory(src)
option(ENABLE_EXAMPLES "Build examples" ON)
if(ENABLE_EXAMPLES)
    add_subdirectory(examples)
    add_subdirectory(userExamples)
endif()

# ---------- Testing (optional) ----------
option(ENABLE_TESTING "Build tests" OFF)
if(ENABLE_TESTING)
    add_subdirectory(tests)
endif()

 # ---------- Python bindings (optional) ----------
option(BUILD_PYBINDINGS "Build Python bindings (pyOILAB)" ON)
if(BUILD_PYBINDINGS)
    add_subdirectory(bindings)
endif()

# ---------- Installation ----------
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# Install header files (preserves subdirectory structure under oILAB/)
install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/oILAB
    FILES_MATCHING PATTERN "*.h"
)

# Generate package config files
configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/oILABConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/oILABConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/oILAB
)

set(OILAB_VERSION_NUMBER "${PROJECT_VERSION}")
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/oILABConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

# Install CMake config files
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/oILABConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/oILABConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/oILAB
)

# Export targets
install(EXPORT oILABTargets
    FILE oILABTargets.cmake
    NAMESPACE oILAB::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/oILAB
)

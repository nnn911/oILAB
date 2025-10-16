add_library(oILAB SHARED    Lattices/BiCrystal.cpp
                            Lattices/Gb.cpp
                            Lattices/Lattice.cpp 
                            Lattices/LatticeVector.cpp
                            Math/RationalMatrix.cpp
                            Lattices/LatticeCore.cpp
                            Math/RLLL.cpp
                            Lattices/ReciprocalLatticeDirection.cpp
                            Math/LLL.cpp
                            Lattices/LatticeDirection.cpp
                            Lattices/RationalLatticeDirection.cpp
                            Lattices/RationalReciprocalLatticeDirection.cpp
                            Lattices/ReciprocalLatticeVector.cpp
                            Math/Farey.cpp
                            Lattices/GbShifts.cpp
                            Lattices/GbMaterialTensors.cpp)


# Conditionally apply the export property for MSVC on Windows
if(MSVC)
    set_property(TARGET oILAB PROPERTY WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()

target_include_directories(oILAB
    PUBLIC
        # For build tree: point to include/ and its subdirectories for current includes to work
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/IO>
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/Lattices>
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/MonteCarlo>
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/Math>
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/Utilities>
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/Visualization>
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/Bindings>
        # For install tree: point to oILAB/ so users can do #include <oILAB/Lattices/LatticeModule.h>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
        # Also add subdirectories for compatibility with non-prefixed includes
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/oILAB/IO>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/oILAB/Lattices>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/oILAB/MonteCarlo>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/oILAB/Math>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/oILAB/Utilities>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/oILAB/Visualization>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/oILAB/Bindings>
)

target_link_libraries(oILAB PUBLIC Eigen3::Eigen)

# Install library
install(TARGETS oILAB
    EXPORT oILABTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)


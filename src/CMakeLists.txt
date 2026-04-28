add_library(oILAB INTERFACE)

target_include_directories(oILAB
    INTERFACE
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

target_link_libraries(oILAB INTERFACE Eigen3::Eigen)

install(TARGETS oILAB
    EXPORT oILABTargets
)


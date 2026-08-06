# add the executable
add_executable(testGB3dMixed testGB3dMixed.cpp)
target_link_libraries(testGB3dMixed oILAB)

add_test(TestGB3dMixed testGB3dMixed)

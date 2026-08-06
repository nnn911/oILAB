# add the executable
add_executable(testGB3dTwist testGB3dTwist.cpp)
target_link_libraries(testGB3dTwist oILAB)

add_test(TestGB3dTwist testGB3dTwist)

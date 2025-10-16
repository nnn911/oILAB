# oILAB

Documentation
-------------
https://oilab-project.github.io/oILAB/


Requirements: Eigen 3.4 or later versions

To Build
--------

1) mkdir build
2) cd build
3) cmake -DCMAKE_BUILD_TYPE=Release ..
4) make (or make -j 8 for faster parallel compilation)

To Install
----------

After building, you can install oILAB system-wide or to a custom location:

### System-wide installation (requires sudo):
```bash
cd build
sudo make install
```

### Custom installation location:
```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/path/to/install ..
make -j 8
make install
```

### Using oILAB in your CMake project:
After installation, you can use oILAB in your CMake project by adding:
```cmake
find_package(oILAB REQUIRED)
target_link_libraries(your_target PRIVATE oILAB::oILAB)
```

In your C++ code, include headers using:
```cpp
#include <oILAB/Lattices/LatticeModule.h>
#include <oILAB/Math/RationalMatrix.h>
// etc.
```

If you installed to a custom location, set the installation path:
```bash
cmake -DCMAKE_PREFIX_PATH=/path/to/install ..
```


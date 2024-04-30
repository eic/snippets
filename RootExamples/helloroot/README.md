### Helloroot

A small sekeleton to demonstrate how to link standalone C++ code to ROOT libraries. To use, rename and adapt the source files and `CMakeLists.txt` to your project.

To build using cmake, create a build directory, navigate to it and run cmake. e.g.:

```
mkdir build
cd build
cmake .. 
make 
```
You can specify a number of parallel build threads with the -j flag, e.g.
```
make -j4
```

You can specify an install directory to cmake with
-DCMAKE_INSTALL_PREFIX=<path>
then, after building, 
```
make install
```
to install the headers and libraries under that location.
There is no "make uninstall" but (on Unix-like systems)
you can do
```
xargs rm < install_manifest.txt
```
from the cmake build directory.

Note: The particular cmake configuration is designed to automatically use all the compile options needed for root, e.g. the c++ std version. On some systems this mechanism fails and you may have to uncomment l. 30 or l. 31 in `CMakeLists.txt`.



# CMakeLists.txt for helloroot.
# More complicated than needed but demonstrates making and linking your own libraries
# cf. https://cliutils.gitlab.io/modern-cmake/chapters/packages/ROOT.html
# https://root.cern/manual/integrate_root_into_my_cmake_project/

cmake_minimum_required(VERSION 3.10)
project(helloroot VERSION 1.0 LANGUAGES CXX ) # not needed

 # Find ROOT. Use at least 6.20 for smoother cmake support
 find_package(ROOT 6.20 REQUIRED )

 message ( " ROOT Libraries = " ${ROOT_LIBRARIES} )

 ###############################################################################################

# Main target is the libhelloroot library
add_library(
  # You can use wildcards but it's cleaner to list the files explicitly
   helloroot
   SHARED
   src/helloroot.cxx
   )

# The particular syntax here is a bit annoying because you have to list all the sub-modules you need
# but it picks up automatically all the compile options needed for root, e.g. the c++ std version
# You can find all available ROOT modules with `root-config --libs`
target_link_libraries(helloroot PUBLIC ROOT::Core ROOT::RIO ROOT::Rint ROOT::Tree ROOT::EG ROOT::Physics )

## The above _should_ be true, and it is on most systems. If it's not, uncoment one of the following lines
# target_compile_features(helloroot PUBLIC cxx_std_17)
# target_compile_features(helloroot PUBLIC cxx_std_20)

# include directories - this is also overkill but useful if you want to create dictionaries
# Contact kkauder@gmail.com for that - it's too much for this example
target_include_directories(helloroot 
PUBLIC 
$<INSTALL_INTERFACE:include>
$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
 )

# Can add addtional options here
target_compile_options(helloroot PRIVATE -Wall -Wextra -pedantic -g)  

##############################################################################################################

## Build executables
add_executable(helloexec src/helloexec.cxx)
# target_compile_options(helloexec PRIVATE -Wall -Wextra -pedantic -g)
target_link_libraries(helloexec helloroot )
target_include_directories(helloexec
  PRIVATE
  ${ROOT_INCLUDE_DIRS}
  )

install(TARGETS helloexec DESTINATION bin)


##############################################################################################################

## Install library
# Could also use include(GNUInstallDirs)
# and then destinations of the form ${CMAKE_INSTALL_INCLUDEDIR}
install(TARGETS helloroot
  EXPORT helloroot-export
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  )

## Install headers
install (DIRECTORY ${CMAKE_SOURCE_DIR}/include/helloroot
  DESTINATION  include/helloroot
  )

## Generate configuration file - this allows you to use cmake in another project 
## to find and link the installed helloroot library
install(EXPORT helloroot-export
  FILE
  hellorootConfig.cmake
  NAMESPACE
    helloroot::
  DESTINATION
  cmake
  )

## Final message
message( " Done!")



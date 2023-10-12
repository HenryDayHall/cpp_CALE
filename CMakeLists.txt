cmake_minimum_required(VERSION 3.4...3.18)
project(Sgwj_library)

# Useful definitions
set(PROJECT_HEADER_DIR ${CMAKE_SOURCE_DIR}/sgwj)
set(PROJECT_SOURCE_DIR ${CMAKE_SOURCE_DIR}/src)

set(PROJECT_HEADERS 
    ${PROJECT_HEADER_DIR}/cluster.hxx
    ${PROJECT_HEADER_DIR}/functions.hxx
   )
set(PROJECT_SOURCES
    ${PROJECT_SOURCE_DIR}/cluster.cxx
    ${PROJECT_SOURCE_DIR}/functions.cxx
   )
message(PROJECT_NAME: ${PROJECT_NAME})
message(PROJECT_HEADERS: ${PROJECT_HEADERS})
message(PROJECT_SOURCES: ${PROJECT_SOURCES})

# Add a library target
add_library(${PROJECT_NAME} STATIC ${PROJECT_HEADERS} ${PROJECT_SOURCES})

# Include directories for the library
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_HEADER_DIR})

# Add an executable target
add_executable(Sgwj_toy examples/toy.cxx)

# Link it to the library
target_link_libraries(Sgwj_toy PRIVATE Sgwj_library)

# Now create python bindings
project(Sgwj_pybind)
add_subdirectory(pybind11)
pybind11_add_module(Sgwj_pybind src/bindings.cxx)
target_link_libraries(Sgwj_pybind PRIVATE Sgwj_library)
target_compile_definitions(Sgwj_pybind
                           PRIVATE VERSION_INFO=0.0.1)

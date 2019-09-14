cmake_minimum_required (VERSION 2.8.11)
set (CMAKE_CXX_STANDARD 11)

project (sine)

set (CMAKE_CXX_FLAGS "-Wno-deprecated ${CMAKE_CXX_FLAGS}")

add_definitions(-DBUILD_COMPONENT_SRC_PREFIX="" -DBUILD_OPTLEVEL_DEV)

set(USD_ROOT_DIR ~/USD)
set(USD_LIB_DIR ${USD_ROOT_DIR}/lib/)
set(USD_INC_DIR ${USD_ROOT_DIR}/include/)

add_executable (sine
        main.cpp
        )

find_package(PythonLibs)

target_include_directories(sine
        PUBLIC
            ${USD_INC_DIR}
            ${PYTHON_INCLUDE_PATH}
        )

target_link_libraries(
        sine
                ${USD_LIB_DIR}/libgf.so
                ${USD_LIB_DIR}/libsdf.so
                ${USD_LIB_DIR}/libtf.so
                ${USD_LIB_DIR}/libusd.so
                ${USD_LIB_DIR}/libusdGeom.so
                ${USD_LIB_DIR}/libvt.so
                ${USD_LIB_DIR}/libboost_python.so # todo find library
                ${PYTHON_LIBRARY}
        )
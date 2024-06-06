cmake_minimum_required(VERSION 3.27)
project(OpenMP)

enable_language(CXX)

if (MSVC)
    # Замените компилятор на g++
    set(CMAKE_CXX_COMPILER g++)
    set(CMAKE_CXX_STANDARD 17)

    find_package(OpenMP REQUIRED)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -fopenmp")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
else()
    # Если не MSVC, используйте g++
    set(CMAKE_CXX_COMPILER g++)
    set(CMAKE_CXX_STANDARD 17)

    find_package(OpenMP REQUIRED)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

set(SOURCE_FILES_1 1.cpp)
add_executable(1 ${SOURCE_FILES_1})

set(SOURCE_FILES_2 2.cpp)
add_executable(2 ${SOURCE_FILES_2})

set(SOURCE_FILES_3 3.cpp)
add_executable(3 ${SOURCE_FILES_3})

set(SOURCE_FILES_4 4.cpp)
add_executable(4 ${SOURCE_FILES_4})

set(SOURCE_FILES_5 5.cpp)
add_executable(5 ${SOURCE_FILES_5})

set(SOURCE_FILES_6 6.cpp)
add_executable(6 ${SOURCE_FILES_6})

set(SOURCE_FILES_7 7.cpp)
add_executable(7 ${SOURCE_FILES_7})

set(SOURCE_FILES_10 10.cpp)
add_executable(10 ${SOURCE_FILES_10})

set(SOURCE_FILES_11 11.cpp)
add_executable(11 ${SOURCE_FILES_11})
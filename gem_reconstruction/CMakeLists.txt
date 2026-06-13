cmake_minimum_required(VERSION 3.16)
project(gem LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# target
add_library(gem SHARED
    src/GEMCluster.cpp
    src/GEMPlane.cpp
    src/Cuts.cpp
    src/ValueType.cpp
    src/GEMDetector.cpp
    src/GEMDetectorLayer.cpp
    src/MPDDataStruct.cpp
    third_party/ConfigArgs.cpp
    third_party/ConfigObject.cpp
    third_party/ConfigOption.cpp
    third_party/ConfigParser.cpp
    third_party/ConfigValue.cpp
)

# Include directories
target_include_directories(gem
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/third_party
)


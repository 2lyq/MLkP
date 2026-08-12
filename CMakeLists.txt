cmake_minimum_required(VERSION 4.0)
project(Regional_MSRs)

set(CMAKE_CXX_STANDARD 20)

add_executable(Regional_MSRs main.cpp Types.h FindNeibs.h FindNeibs.cpp Model2.h Model2.cpp FindRegion.h FindRegion.cpp CSVReader/CSVReader.h CSVReader/CSVReader.cpp EventClock/EventClock.hpp)

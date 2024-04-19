cmake_minimum_required(VERSION 3.1...3.28)

project(
  seamplusplus
  VERSION 1.0
  LANGUAGES CXX)

file(DOWNLOAD https://raw.githubusercontent.com/lvandeve/lodepng/master/lodepng.cpp src/lodepng.cpp)

file(DOWNLOAD https://raw.githubusercontent.com/lvandeve/lodepng/master/lodepng.h src/lodepng.h)

add_executable(seamplusplus src/main.cpp src/lodepng.cpp)

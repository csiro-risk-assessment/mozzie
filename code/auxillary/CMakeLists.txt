cmake_minimum_required(VERSION 3.6)
project(ab_convert)
add_executable(ab_convert ab_convert.c ../csvparser.c)
target_compile_options(ab_convert PUBLIC -O3 -fno-strict-aliasing)

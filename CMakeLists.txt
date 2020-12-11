cmake_minimum_required(VERSION 3.9)
project(twistedHessian)

include(FindPkgConfig)

pkg_check_modules(GMP REQUIRED IMPORTED_TARGET gmp)

add_executable(${PROJECT_NAME} "main.c" "funcs.h" "funcs.c")

#find_package(GMP REQUIRED)

target_link_libraries(${PROJECT_NAME} PRIVATE PkgConfig::GMP)





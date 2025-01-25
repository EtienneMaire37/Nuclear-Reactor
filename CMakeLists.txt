cmake_minimum_required(VERSION 3.10)

# Set the project name
project(NuclearReactor)

# Specify the C standard
set(CMAKE_C_STANDARD 99)

# Set the C compiler to gcc
set(CMAKE_C_COMPILER gcc)

# Set the C flags
set(CMAKE_C_FLAGS "-Wall -Wextra -Werror -Wpedantic -Ofast -lm")

# Add the executable
add_executable(NuclearReactor src/main.c)
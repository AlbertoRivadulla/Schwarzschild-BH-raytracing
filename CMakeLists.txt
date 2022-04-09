cmake_minimum_required(VERSION 3.7)

# Obtain a file compile_commands.json used by ccls (through the plugin coc.nvim)
# to provide code completion in Neovim
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE INTERNAL "")

# Name of the project
project(project)

# set(LIB_DIR "${PROJECT_SOURCE_DIR}/lib")

# Find the SDL2 library
# find_package(SDL2 REQUIRED)

# Create a variable with all the include directories
# Add SDL2 into include directories
set(INCLUDE
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/src
    # ${PROJECT_SOURCE_DIR}/include
    # ${SDL2_INCLUDE_DIRS}
   )

include_directories(${INCLUDE})

# Create a variable with a link to all cpp files to compile
set(SOURCES
    ${PROJECT_SOURCE_DIR}/main.cpp
    ${PROJECT_SOURCE_DIR}/src/rungeKutta.cpp
    ${PROJECT_SOURCE_DIR}/src/saveToFiles.cpp
    ${PROJECT_SOURCE_DIR}/src/traceRays.cpp
   )

add_executable(main ${SOURCES})

# target_link_libraries(main SDL2 SDL2_image SDL2_ttf SDL2_mixer)

# Get rid of the cmake_install.cmake file created
set(CMAKE_SKIP_INSTALL_RULES True)

# # STB
# set(STB_DIR "${LIB_DIR}/stb")
# target_include_directories(main PRIVATE "${STB_DIR}")

# Move the file compile_commands.json to the directory of the project.
execute_process(COMMAND cp compile_commands.json ..) 

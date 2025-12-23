cmake_minimum_required(VERSION 3.10)
project(FlagsK4FreePartitions CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

set(COMMON_FLAG_SOURCES
  src/option-parser.cpp
  src/solution-processor.cpp
  src/flag-generator.cpp
  src/flag.cpp
  src/utils.cpp
  src/cache.cpp
  src/solver.cpp
)

# 2-partition, simple proof for the case when minimum degree is at least 1/2.
add_executable(sum-2 problems/sum-2.cpp ${COMMON_FLAG_SOURCES})
target_include_directories(sum-2 PRIVATE src)
target_compile_definitions(sum-2 PRIVATE V=6 G_COLORED_EDGES=2)

set(sum-2-sdp-directory sum-2)

add_custom_target(sum-2-solve-approximate
  COMMAND ./bin/sum-2 -l 4 -s ${sum-2-sdp-directory}
  DEPENDS sum-2
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

add_custom_target(sum-2-round
  COMMAND sage round.sage --sdp-directory ${sum-2-sdp-directory} --scaling-factor 540
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

# 2-partition, general case proof leveraging colored vertices (hence color in the names :)).
add_executable(sum-2-color problems/sum-2-color.cpp ${COMMON_FLAG_SOURCES})
target_include_directories(sum-2-color PRIVATE src)
target_compile_definitions(sum-2-color PRIVATE V=6 G_COLORED_EDGES=2 G_COLORED_VERTICES=2)

set(sum-2-color-sdp-directory sum-2-color)

add_custom_target(sum-2-color-solve-approximate
  COMMAND ./bin/sum-2-color -l 4 -s ${sum-2-color-sdp-directory}
  DEPENDS sum-2-color
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

add_custom_target(sum-2-color-round
  COMMAND sage round.sage --sdp-directory ${sum-2-color-sdp-directory} --scaling-factor 540
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

add_executable(max-2 problems/max-2.cpp ${COMMON_FLAG_SOURCES})
target_include_directories(max-2 PRIVATE src)
target_compile_definitions(max-2 PRIVATE V=5 G_COLORED_EDGES=2 G_COLORED_VERTICES=4)

add_executable(max-2-less-colors problems/max-2.cpp ${COMMON_FLAG_SOURCES})
target_include_directories(max-2-less-colors PRIVATE src)
target_compile_definitions(max-2-less-colors PRIVATE V=6 G_COLORED_EDGES=2 G_COLORED_VERTICES=3)

set(max-2-case-1-sdp-directory max-2-case-1)
add_custom_target(max-2-case-1-solve-approximate
  COMMAND ./bin/max-2 -l 3 -s ${max-2-case-1-sdp-directory} -c 1
  DEPENDS max-2
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

set(max-2-case-2-sdp-directory max-2-case-2)
add_custom_target(max-2-case-2-solve-approximate
  COMMAND ./bin/max-2 -l 3 -s ${max-2-case-2-sdp-directory} -c 2
  DEPENDS max-2
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

add_executable(max-2-experiment problems/max-2-experiment.cpp ${COMMON_FLAG_SOURCES})
target_include_directories(max-2-experiment PRIVATE src)
target_compile_definitions(max-2-experiment PRIVATE V=6 G_COLORED_EDGES=2 G_COLORED_VERTICES=2)

set(max-2-experiment-sdp-directory max-2-experiment)

add_custom_target(max-2-experiment-solve-approximate
  COMMAND ./bin/max-2-experiment -l 4 -s ${max-2-experiment-sdp-directory}
  DEPENDS max-2-experiment
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

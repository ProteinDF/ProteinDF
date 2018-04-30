include_directories(${CMAKE_CURRENT_SOURCE_DIR}/libpdftl)

add_subdirectory(libpdftl)
add_subdirectory(libpdf)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/libpdf)
add_subdirectory(pdf)

add_subdirectory(tools)


if (GTEST_FOUND)
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/unit_test)
  add_subdirectory(unit_test)
endif (GTEST_FOUND)

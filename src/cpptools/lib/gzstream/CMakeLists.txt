add_library(gzstream STATIC gzstream.C)
target_include_directories(gzstream PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(gzstream z)

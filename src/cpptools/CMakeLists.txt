cmake_minimum_required(VERSION 3.10)
project(MyTools)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-O3")

# === Check for Conda environment ===
if(NOT DEFINED ENV{CONDA_PREFIX})
    message(FATAL_ERROR "Conda environment not activated! Please run:\n  conda activate pangaea")
endif()

message(STATUS "Using Conda environment at: $ENV{CONDA_PREFIX}")

# === Set Conda-specific paths ===
set(CONDA_PREFIX "$ENV{CONDA_PREFIX}")
set(CONDA_LIB_DIR "${CONDA_PREFIX}/lib")
set(CONDA_INCLUDE_DIR "${CONDA_PREFIX}/include")

# Make sure CMake searches Conda directories
list(APPEND CMAKE_PREFIX_PATH "${CONDA_PREFIX}")
list(APPEND CMAKE_LIBRARY_PATH "${CONDA_LIB_DIR}")
list(APPEND CMAKE_INCLUDE_PATH "${CONDA_INCLUDE_DIR}")

# Make linker prefer Conda libs
link_directories("${CONDA_LIB_DIR}")

# === Find required packages ===
find_package(OpenMP REQUIRED)

find_library(HTS_LIB hts PATHS "${CONDA_LIB_DIR}" NO_DEFAULT_PATH REQUIRED)
find_library(ZLIB_LIB z PATHS "${CONDA_LIB_DIR}" NO_DEFAULT_PATH REQUIRED)
find_library(BZ2_LIB bz2 PATHS "${CONDA_LIB_DIR}" NO_DEFAULT_PATH REQUIRED)

find_path(HTS_INCLUDE_DIR hts.h PATHS "${CONDA_INCLUDE_DIR}/htslib" NO_DEFAULT_PATH REQUIRED)
find_path(ZLIB_INCLUDE_DIR zlib.h PATHS "${CONDA_INCLUDE_DIR}" NO_DEFAULT_PATH REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${ZLIB_INCLUDE_DIR}")

# === Output directory for binaries ===
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/../bin")

# === Add gzstream library ===
add_subdirectory(lib/gzstream)

# === Common include paths ===
set(COMMON_INCLUDES
    "${CMAKE_SOURCE_DIR}/lib/cmdline"
    "${CMAKE_SOURCE_DIR}/lib/gzstream"
    "${CMAKE_SOURCE_DIR}/lib/ThreadPool"
    "${HTS_INCLUDE_DIR}"
    "${ZLIB_INCLUDE_DIR}"
)

# === Common libraries ===
set(COMMON_LIBS pthread gzstream)

# === Define executables ===

# count_kmer
add_executable(count_kmer count_kmer.cpp)
target_include_directories(count_kmer PRIVATE ${COMMON_INCLUDES})
target_link_libraries(count_kmer PRIVATE ${COMMON_LIBS})

# count_tnf
add_executable(count_tnf count_tnf.cpp)
target_include_directories(count_tnf PRIVATE ${COMMON_INCLUDES})
target_link_libraries(count_tnf PRIVATE ${COMMON_LIBS})

# extract_reads
add_executable(extract_reads extract_reads.cpp)
target_include_directories(extract_reads PRIVATE ${COMMON_INCLUDES})
target_link_libraries(extract_reads PRIVATE ${COMMON_LIBS})

# parse_header
add_executable(parse_header parse_header.cpp)
target_include_directories(parse_header PRIVATE ${COMMON_INCLUDES})
target_link_libraries(parse_header PRIVATE ${COMMON_LIBS})

# assign_barcodes
add_executable(assign_barcodes assign_barcodes.cpp)
target_include_directories(assign_barcodes PRIVATE ${COMMON_INCLUDES})
target_link_libraries(assign_barcodes PRIVATE ${COMMON_LIBS})
target_compile_options(assign_barcodes PRIVATE -fopenmp)
target_link_options(assign_barcodes PRIVATE -fopenmp)

# extract_unmapped
add_executable(extract_unmapped extract_unmapped.cpp)
target_include_directories(extract_unmapped PRIVATE ${COMMON_INCLUDES})
target_link_libraries(extract_unmapped PRIVATE ${COMMON_LIBS} ${HTS_LIB} ${BZ2_LIB})
target_compile_options(extract_unmapped PRIVATE -fopenmp)
target_link_options(extract_unmapped PRIVATE -fopenmp)

# add_barcode
add_executable(add_barcode add_barcode.cpp)
target_include_directories(add_barcode PRIVATE ${COMMON_INCLUDES})
target_link_libraries(add_barcode PRIVATE ${COMMON_LIBS} ${HTS_LIB})
target_compile_options(add_barcode PRIVATE -fopenmp)
target_link_options(add_barcode PRIVATE -fopenmp)

# preprocess_stlfr
add_executable(preprocess_stlfr preprocess_stlfr.cpp)
target_include_directories(preprocess_stlfr PRIVATE ${COMMON_INCLUDES})
target_link_libraries(preprocess_stlfr PRIVATE ${COMMON_LIBS})

# preprocess_tellseq
add_executable(preprocess_tellseq preprocess_tellseq.cpp)
target_include_directories(preprocess_tellseq PRIVATE ${COMMON_INCLUDES})
target_link_libraries(preprocess_tellseq PRIVATE ${COMMON_LIBS})
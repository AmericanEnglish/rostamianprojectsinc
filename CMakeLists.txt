cmake_minimum_required(VERSION 3.10)

project(RostamianProjects)
enable_language(C)
set("CMAKE_C_COMPILER" "gcc")
set("XMALLOC" "${PROJECT_SOURCE_DIR}/xmalloc/")
set("ARRAY" "${PROJECT_SOURCE_DIR}/array/")
set("XMALLOCSRC" "${XMALLOC}/xmalloc.c")
#target_include_directories("xmalloc/", "array/")
include_directories("${PROJECT_SOURCE_DIR}/xmalloc/" "${PROJECT_SOURCE_DIR}/array/")


add_executable("xmalloc-demo-1" "${XMALLOC}/xmalloc-demo-1.c" "${XMALLOCSRC}")
add_executable("xmalloc-demo-2" "${XMALLOC}/xmalloc-demo-2.c" "${XMALLOCSRC}")
add_executable("hilbert_matrix" "${ARRAY}/hilbert_matrix.c" ${XMALLOCSRC})


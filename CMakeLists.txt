cmake_minimum_required(VERSION 3.10)

project(RostamianProjects)
enable_language(C)
set("CMAKE_C_COMPILER" "gcc")

# Source directories
set("XMALLOC" "${PROJECT_SOURCE_DIR}/xmalloc/")
set("XMALLOCSRC" "${XMALLOC}/xmalloc.c")

set("ARRAY" "${PROJECT_SOURCE_DIR}/array/")

set("VERTS" "${PROJECT_SOURCE_DIR}/vertices/")
set("VERTSRC" "${VERTS}/rank-vertices.c")

set("NMEAD" "${PROJECT_SOURCE_DIR}/nelder-mead/")
set("NMEADSRC" "${NMEAD}/nelder-mead.c")
# Header directories
include_directories("${XMALLOC}" "${ARRAY}" "${VERTS}" "${NMEAD}")
set("DEMO2DSRC" "${NMEADSRC}" "${XMALLOCSRC}")

set("CMAKE_C_FLAGS" "-Wall -pedantic -std=c99 -O2")
add_executable("xmalloc-demo-1" "${XMALLOC}/xmalloc-demo-1.c" "${XMALLOCSRC}")
add_executable("xmalloc-demo-2" "${XMALLOC}/xmalloc-demo-2.c" "${XMALLOCSRC}")
add_executable("hilbert_matrix" "${ARRAY}/hilbert_matrix.c" ${XMALLOCSRC})
add_executable("rank-vertices-demo" "${VERTS}/rank-vertices-demo.c" )
add_executable("demo2D" "${NMEAD}/demo2D.c" "${DEMO2DSRC}")

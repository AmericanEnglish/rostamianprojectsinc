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

set("NN" "${PROJECT_SOURCE_DIR}/neuralnet/")
set("NNSRC" "${NN}/neural-net-ode.c")

set("DEMO2DSRC" "${NMEADSRC}" "${XMALLOCSRC}")
set("DEMOENSRC" "${NMEADSRC}" "${XMALLOCSRC}")
set("DEMOODESRC" "${NNSRC}" "${XMALLOCSRC}" "${NMEADSRC}")

# Header directories
include_directories("${XMALLOC}" "${ARRAY}" "${VERTS}" "${NMEAD}" "${NN}")

set("CMAKE_C_FLAGS" "-Wall -pedantic -std=c99 -O2")
add_executable("xmalloc-demo-1" "${XMALLOC}/xmalloc-demo-1.c" "${XMALLOCSRC}")
add_executable("xmalloc-demo-2" "${XMALLOC}/xmalloc-demo-2.c" "${XMALLOCSRC}")
add_executable("hilbert_matrix" "${ARRAY}/hilbert_matrix.c" ${XMALLOCSRC})
add_executable("rank-vertices-demo" "${VERTS}/rank-vertices-demo.c" )
add_executable("demo2D" "${NMEAD}/demo2D.c" "${DEMO2DSRC}")
add_executable("demo-energy" "${NMEAD}/demo-energy.c" "${DEMOENSRC}")
add_executable("demo-ode1" "${NN}/demo-ode1.c" "${DEMOODESRC}")
target_link_libraries("demo-energy" "m")
target_link_libraries("demo-ode1" "m")

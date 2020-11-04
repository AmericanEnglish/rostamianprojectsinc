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
set("NNSRC" "${NN}/neural-nets-ode.c")

# The source for these files can be seen here:
# https://userpages.umbc.edu/~rostamia/cbook/neural-nets-odes.html
# We will download these right here for completion.
if (NOT EXISTS "${NN}/plot-with-maple.c")
    file(DOWNLOAD "https://userpages.umbc.edu/~rostamia/cbook/neural-nets-odes/plot-with-maple.c" "${NN}/plot-with-maple.c" 
        SHOW_PROGRESS) 
endif()

if (NOT EXISTS "${NN}/plot-with-matlab.c")
    file(DOWNLOAD "https://userpages.umbc.edu/~rostamia/cbook/neural-nets-odes/plot-with-matlab.c" "${NN}/plot-with-matlab.c" 
        SHOW_PROGRESS) 
endif()
set("ROSTPLOTSRC" "${NN}/plot-with-matlab.c" "${NN}/plot-with-maple.c")

set("DEMO2DSRC" "${NMEADSRC}" "${XMALLOCSRC}")
set("DEMOENSRC" "${NMEADSRC}" "${XMALLOCSRC}")
set("DEMOODESRC" "${NNSRC}" "${XMALLOCSRC}" "${NMEADSRC}" "${ROSTPLOTSRC}")

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

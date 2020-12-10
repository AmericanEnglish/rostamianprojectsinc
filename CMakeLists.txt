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

set("NNO" "${PROJECT_SOURCE_DIR}/neuralnet-ode/")
set("NNOSRC" "${NNO}/neural-nets-ode.c")

set("NNP"    "${PROJECT_SOURCE_DIR}/neuralnet-pde/")
set("NNPSRC" "${PROJECT_SOURCE_DIR}/neuralnet-pde/neural-nets-pde.c")

set("FD"      "${PROJECT_SOURCE_DIR}/finite-differences/")
set("HEATSRC" "${PROJECT_SOURCE_DIR}/finite-differences/heat-solve.c")

# The source for these files can be seen here:
# https://userpages.umbc.edu/~rostamia/cbook/neural-nets-odes.html
# We will download these right here for completion.
if (NOT EXISTS "${NNO}/plot-with-maple.c")
    file(DOWNLOAD "https://userpages.umbc.edu/~rostamia/cbook/neural-nets-odes/plot-with-maple.c" "${NNO}/plot-with-maple.c" SHOW_PROGRESS) 
endif()

if (NOT EXISTS "${NNO}/plot-with-matlab.c")
    file(DOWNLOAD "https://userpages.umbc.edu/~rostamia/cbook/neural-nets-odes/plot-with-matlab.c" "${NNO}/plot-with-matlab.c" SHOW_PROGRESS) 
endif()
set("ROSTPLOTODESRC" "${NNO}/plot-with-matlab.c" "${NNO}/plot-with-maple.c")

# The source for these files can be seen here:
# https://userpages.umbc.edu/~rostamia/cbook/neural-nets-pdes.html
# We will download these right here for completion.
if (NOT EXISTS "${NNP}/plot-with-maple.c")
    file(DOWNLOAD "https://userpages.umbc.edu/~rostamia/cbook/neural-nets-pdes/plot-with-maple.c"
        "${NNP}/plot-with-maple.c" SHOW_PROGRESS) 
endif()

if (NOT EXISTS "${NNP}/plot-with-matlab.c")
    file(DOWNLOAD
        "https://userpages.umbc.edu/~rostamia/cbook/neural-nets-pdes/plot-with-matlab.c"
        "${NNP}/plot-with-matlab.c" SHOW_PROGRESS) 
endif()
set("ROSTPLOTPDESRC" "${NNP}/plot-with-matlab.c" "${NNP}/plot-with-maple.c")

set("DEMO2DSRC" "${NMEADSRC}" "${XMALLOCSRC}")
set("DEMOENSRC" "${NMEADSRC}" "${XMALLOCSRC}")
set("DEMOODESRC" "${NNOSRC}" "${XMALLOCSRC}" "${NMEADSRC}" "${ROSTPLOTODESRC}")

# The PDE program has a lot of demos
set("DEMOSQSRC" "${NNPSRC}" "${XMALLOCSRC}" "${NMEADSRC}" "${ROSTPLOTPDESRC}") 
set("DEMOLSRC" "${NNPSRC}" "${XMALLOCSRC}" "${NMEADSRC}" "${ROSTPLOTPDESRC}") 


# This downloads the 3d plotters for the rostamian finite difference code
if (NOT EXISTS "${FD}/plot3d.c")
    file(DOWNLOAD "https://userpages.umbc.edu/~rostamia/cbook/fd1/plot3d.c" 
        "${FD}/plot3d.c" SHOW_PROGRESS) 
endif()
if (NOT EXISTS "${FD}/plot3d.h")
    file(DOWNLOAD "https://userpages.umbc.edu/~rostamia/cbook/fd1/plot3d.h" 
        "${FD}/plot3d.h" SHOW_PROGRESS) 
endif()
set("ROSTPLOTFDSRC" "${FD}/plot3d.c")

# The finite difference program 
set("FDSRC" "${HEATSRC}" "${XMALLOCSRC}" "${ROSTPLOTFDSRC}")

# Header directories
include_directories("${XMALLOC}" "${ARRAY}" "${VERTS}" "${NMEAD}" "${NNO} ${FD}")


set("CMAKE_C_FLAGS" "-Wall -pedantic -std=c99 -O2")
add_executable("xmalloc-demo-1" "${XMALLOC}/xmalloc-demo-1.c" "${XMALLOCSRC}")
add_executable("xmalloc-demo-2" "${XMALLOC}/xmalloc-demo-2.c" "${XMALLOCSRC}")
add_executable("hilbert_matrix" "${ARRAY}/hilbert_matrix.c" ${XMALLOCSRC})
add_executable("rank-vertices-demo" "${VERTS}/rank-vertices-demo.c" )
add_executable("demo2D" "${NMEAD}/demo2D.c" "${DEMO2DSRC}")
add_executable("demo-energy" "${NMEAD}/demo-energy.c" "${DEMOENSRC}")
add_executable("demo-ode1" "${NNO}/demo-ode1.c" "${DEMOODESRC}")
add_executable("demo-square" "${NNP}/demo-square.c" "${DEMOSQSRC}")
add_executable("demo-L-shaped" "${NNP}/demo-L-shaped-fragment.c" "${DEMOLSRC}")
add_executable("demo-square-hole" "${NNP}/demo-square-with-hole-fragment.c" "${DEMOLSRC}")
add_executable("demo-finite-methods-1" "${FD}/demo1.c" "${FDSRC}")
add_executable("demo-finite-methods-2" "${FD}/demo2.c" "${FDSRC}")
add_executable("demo-finite-methods-3" "${FD}/demo3.c" "${FDSRC}")
add_executable("demo-finite-methods-4" "${FD}/demo4.c" "${FDSRC}")

target_link_libraries("demo-energy" "m")
target_link_libraries("demo-ode1" "m")
target_link_libraries("demo-square" "m")
target_link_libraries("demo-square-hole" "m")
target_link_libraries("demo-L-shaped" "m")
target_link_libraries("demo-finite-methods-1" "m")
target_link_libraries("demo-finite-methods-2" "m")
target_link_libraries("demo-finite-methods-3" "m")
target_link_libraries("demo-finite-methods-4" "m")

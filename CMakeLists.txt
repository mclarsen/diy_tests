if (NOT EXISTS ${PROJECT_SOURCE_DIR}/blt/SetupBLT.cmake)
    message(FATAL_ERROR
        "The BLT submodule is not present. "
        "Run the following two commands in your git repository: \n"
        "    git submodule init \n "
        "    git submodule update" )
endif()

option(ENABLE_MPI         "Build MPI Support"         ON)

include(blt/SetupBLT.cmake)
add_subdirectory(diy)



blt_add_executable(
    NAME        reduce 
    SOURCES     reduce.cpp 
    DEPENDS_ON  diy mpi)

blt_add_executable(
    NAME        all_to_all 
    SOURCES     all_to_all.cpp 
    DEPENDS_ON  diy mpi)

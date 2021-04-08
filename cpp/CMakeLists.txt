if (SUPERNOVA)
    include_directories(${SC_PATH}/external_libraries/nova-tt)
    # actually just boost.atomic
    include_directories(${SC_PATH}/external_libraries/boost)
endif()

include_directories(${CMAKE_SOURCE_DIR}/)
include_directories(scfft_old/)

include_directories(${SC_PATH}/common)
include_directories(${SC_PATH}/plugin_interface)

include_directories(${SC_PATH}/include/plugin_interface)
include_directories(${SC_PATH}/include/common)
include_directories(${SC_PATH}/server/plugins) # for FFT_UGens.h

# old-style directory layout
include_directories(${SC_PATH}/common/Headers/plugin_interface)
include_directories(${SC_PATH}/common/Headers/common)
include_directories(${SC_PATH}/common/Source/plugins) # for FFT_UGens.h

#-------
# macros
#-------

macro(GET_GCC_VERSION VAR)
    if (CMAKE_COMPILER_IS_GNUCC)
        execute_process(
            COMMAND ${CMAKE_CXX_COMPILER} --version OUTPUT_VARIABLE OUTPUT)
        string(REGEX REPLACE
            ".* ([0-9]\\.[0-9]\\.[0-9]) .*" "\\1" ${VAR} ${OUTPUT})
    endif()
endmacro(GET_GCC_VERSION)


macro(FIND_AND_BUILD_PLUGINS DIR)
    file(GLOB PLUGIN_SRC RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "${DIR}/*.cpp")
    foreach(plugin ${PLUGIN_SRC})
        get_filename_component(basename ${plugin} NAME_WE)
        BUILD_PLUGIN(${basename} ${plugin} "" "${DIR}")
    endforeach(plugin)
endmacro()


# Hacky: OSX plugins should install in PLUGIN_DIR
macro(BUILD_PLUGIN PLUGIN_NAME PLUGIN_SOURCES PLUGIN_LIBS PLUGIN_DIR)
    add_library(${PLUGIN_NAME} MODULE ${PLUGIN_SOURCES})
    list(APPEND plugins ${PLUGIN_NAME})
    target_link_libraries(${PLUGIN_NAME} ${PLUGIN_LIBS})

    if(APPLE)
        install(TARGETS ${PLUGIN_NAME}
                DESTINATION ${INSTALL_DESTINATION}/${PLUGIN_DIR})
    endif()

    if(SUPERNOVA)
        add_library(${PLUGIN_NAME}_supernova MODULE ${PLUGIN_SOURCES})
        target_link_libraries(${PLUGIN_NAME}_supernova ${PLUGIN_LIBS})

        list(APPEND supernova_plugins ${PLUGIN_NAME}_supernova)
        set_property(TARGET ${PLUGIN_NAME}_supernova
                     PROPERTY COMPILE_DEFINITIONS SUPERNOVA)
        if(APPLE)
            install(TARGETS ${PLUGIN_NAME}_supernova
                    DESTINATION ${INSTALL_DESTINATION}/${PLUGIN_DIR})
        endif()
    endif()
endmacro()


set(CMAKE_SHARED_MODULE_PREFIX "")
set(PLUGIN_EXTENSION ".so")

if (APPLE OR WIN32)
    set(CMAKE_SHARED_MODULE_SUFFIX ".scx")
    set(PLUGIN_EXTENSION ".scx")
endif()

set(INSTALL_DESTINATION "berlach_plugins")

BUILD_PLUGIN(AmpTrackers
"AmpTrackers.cpp"
"" "plugins")

BUILD_PLUGIN(APFChain
"APFChain.cpp"
"" "plugins")

BUILD_PLUGIN(BerlachFilters
"BerlachFilters.cpp"
"" "plugins")

BUILD_PLUGIN(BitManip4
"BitManip4.cpp"
"" plugins)

BUILD_PLUGIN(BLBufRd
"BLBufRd.cpp"
"" "plugins")

BUILD_PLUGIN(GranularChorus
"GranularChorus.cpp"
"" "plugins")

BUILD_PLUGIN(LowLatConv
"LowLatConv.cpp"
"" "plugins")

BUILD_PLUGIN(PCMRead
"PCMRead.cpp"
"" "plugins")

BUILD_PLUGIN(PhasorLP
"PhasorLP.cpp"
"" "plugins")

BUILD_PLUGIN(ReverbUGens
"ReverbUGens.cpp"
"" "plugins")

BUILD_PLUGIN(SplineHysterisis
"SplineHysterisis.cpp"
"" "plugins")

BUILD_PLUGIN(WarpIIRLattice
"WarpIIRLattice.cpp"
"" "plugins")

BUILD_PLUGIN(WTO
"WTO.cpp"
"" "plugins")

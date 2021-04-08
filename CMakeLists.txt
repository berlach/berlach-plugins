cmake_minimum_required (VERSION 2.8)

project (berlach-plugins)


set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules
                      ${CMAKE_MODULE_PATH})

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/cmake_uninstall.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
  IMMEDIATE @ONLY)

add_custom_target(uninstall
  "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake")

find_package(SuperCollider3)
if (NOT SC_FOUND)
	message(SEND_ERROR "cannot find SuperCollider3 headers. Set the variable SC_PATH.")
else()
	message(STATUS "Using SC source located at ${SC_PATH}")
endif()

set(NOVA_TT_MISSING_ERROR "The nova-tt source code is missing in \
${SC_PATH}/external_libraries/nova-tt.\n Make sure to point to a valid version \
of SuperCollider's source code (with the help of the SC_PATH variable).\n It's \
easiest to use a release tarball of SuperCollider.\n If you use a git clone, \
make sure you have all submodules by running `git submodule update --init \
--recursive` in the root of SuperCollider's source code repository.")

if (NOVA_DISK_IO)
  if (NOT EXISTS "${SC_PATH}/external_libraries/nova-tt/CMakeLists.txt")
    message(FATAL_ERROR "${NOVA_TT_MISSING_ERROR}")
  endif()
endif()

include("${SC_PATH}/SCVersion.txt")
set(PROJECT_VERSION "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}${PROJECT_VERSION_PATCH}")
message(STATUS "Building plugins for SuperCollider version: ${PROJECT_VERSION}")

option(SUPERNOVA "Build plugins for supernova" OFF)

if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
	set(CMAKE_COMPILER_IS_CLANG 1)
endif()

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_CLANG)
    add_definitions(-fvisibility=hidden)

    include (CheckCCompilerFlag)
    include (CheckCXXCompilerFlag)

    CHECK_C_COMPILER_FLAG(-msse HAS_SSE)
    CHECK_CXX_COMPILER_FLAG(-msse HAS_CXX_SSE)

    if (HAS_SSE)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse")
    endif()
    if (HAS_CXX_SSE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse")
    endif()

    CHECK_C_COMPILER_FLAG(-msse2 HAS_SSE2)
    CHECK_CXX_COMPILER_FLAG(-msse2 HAS_CXX_SSE2)

    if (HAS_SSE2)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse2")
    endif()
    if (HAS_CXX_SSE2)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse2")
    endif()

    CHECK_C_COMPILER_FLAG(-mfpmath=sse HAS_FPMATH_SSE)
    CHECK_CXX_COMPILER_FLAG(-mfpmath=sse HAS_CXX_FPMATH_SSE)

    if (HAS_FPMATH_SSE)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mfpmath=sse")
    endif()
    if (HAS_CXX_FPMATH_SSE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mfpmath=sse")
    endif()

	if(NATIVE)
		add_definitions(-march=native)
	endif()

	if(CPP11)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
		if(CMAKE_COMPILER_IS_CLANG)
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
		endif()
	endif()
endif()
if(MINGW)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mstackrealign")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mstackrealign")
endif()


# osx `make install' defaults into cmake_build/SC3-plugins directory
if (APPLE AND IN_PLACE_BUILD)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}")
endif()

# when packaging OSX dmg, install in 'cmake_build/build_osx' so we can include
# directly the created SC3plugins subfolder in the root of the dmg (instead of
# all plugins directories in the dmg root)
if (APPLE AND OSX_PACKAGE)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/build_osx")
endif()

add_subdirectory(cpp)

if (QUARKS OR OSX_PACKAGE) # add quarks DIRECTORY in osx dmg
  add_subdirectory(quarks)
endif()

if (OSX_PACKAGE)
  add_subdirectory(osx_package)
endif()

# Configure readme file
if (WIN32 OR APPLE)
    if (WIN32)
        set(README_FILENAME README_WINDOWS.txt)
    elseif(APPLE)
        set(README_FILENAME README_MACOS.txt)
    else()
        unset(README_FILENAME)
    endif()
endif()

if (DEFINED README_FILENAME)
    configure_file( "${README_FILENAME}.in" "${README_FILENAME}" )
    install( FILES "${CMAKE_CURRENT_BINARY_DIR}/${README_FILENAME}" DESTINATION . RENAME README.txt )
endif()


#############################################
# CPack support

set(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})
set(CPACK_PACKAGE_VERSION ${PROJECT_VERSION})

if(WIN32)
    set(CPACK_GENERATOR ZIP)
    set(CPACK_INCLUDE_TOPLEVEL_DIRECTORY 0)
    set(CPACK_MONOLITHIC_INSTALL 1)
endif()

include(CPack)


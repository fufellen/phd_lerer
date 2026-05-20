cmake_minimum_required(VERSION 3.16)

project(phd_lerer LANGUAGES C)

set(CMAKE_C_STANDARD 90)
set(CMAKE_C_STANDARD_REQUIRED OFF)
set(CMAKE_C_EXTENSIONS ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT WIN32)
    message(FATAL_ERROR "This legacy project currently targets Windows/MSVC because it uses conio.h, kbhit() and getch().")
endif()

if(NOT MSVC)
    message(FATAL_ERROR "Use the Visual Studio/MSVC C compiler. MinGW/Clang may require a compatibility layer for conio.h and old heap APIs.")
endif()

set(POLYNOM_DIR "${CMAKE_CURRENT_SOURCE_DIR}/Полиномы 7_1")
set(NEWLIB_DIR "${CMAKE_CURRENT_SOURCE_DIR}/newLIb")

set(POLYNOM_SOURCES
    "${POLYNOM_DIR}/VIBR.c"
    "${POLYNOM_DIR}/Vvod.c"
    "${POLYNOM_DIR}/COMPL FUN.c"
    "${POLYNOM_DIR}/E_ext.c"
    "${POLYNOM_DIR}/F_df_2.c"
    "${POLYNOM_DIR}/F_df_2_Re.c"
    "${POLYNOM_DIR}/F_df_dop.c"
    "${POLYNOM_DIR}/F_df_IE_MK.c"
    "${POLYNOM_DIR}/FGrin_dop.c"
    "${POLYNOM_DIR}/Fun.c"
    "${POLYNOM_DIR}/Fun mod Bessel.c"
    "${POLYNOM_DIR}/V_.c"
)

set(NEWLIB_SOURCES
    "${NEWLIB_DIR}/CELINEG.C"
    "${NEWLIB_DIR}/FSKCOM.C"
    "${NEWLIB_DIR}/MAT1D.C"
    "${NEWLIB_DIR}/MAT2D.C"
    "${NEWLIB_DIR}/MEMCHK.C"
    "${NEWLIB_DIR}/MUL_DIV.C"
    "${NEWLIB_DIR}/TIM.C"
)

set_source_files_properties(
    ${POLYNOM_SOURCES}
    ${NEWLIB_SOURCES}
    PROPERTIES LANGUAGE C
)

add_executable(phd_lerer
    ${POLYNOM_SOURCES}
    ${NEWLIB_SOURCES}
)

target_include_directories(phd_lerer PRIVATE
    "${POLYNOM_DIR}"
    "${NEWLIB_DIR}"
)

target_compile_definitions(phd_lerer PRIVATE
    _CRT_SECURE_NO_WARNINGS
)

target_compile_options(phd_lerer PRIVATE
    /W3
    /wd4013
    /wd4244
    /wd4305
    /wd4996
)

set_target_properties(phd_lerer PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    RUNTIME_OUTPUT_DIRECTORY_DEBUG "${CMAKE_BINARY_DIR}/bin/Debug"
    RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_BINARY_DIR}/bin/Release"
    RUNTIME_OUTPUT_DIRECTORY_RELWITHDEBINFO "${CMAKE_BINARY_DIR}/bin/RelWithDebInfo"
    RUNTIME_OUTPUT_DIRECTORY_MINSIZEREL "${CMAKE_BINARY_DIR}/bin/MinSizeRel"
)

# R_T_dop.c is not part of the active calculation path and currently duplicates
# globals from F_df_dop.c. Keep it out of this target unless that code path is
# restored deliberately.

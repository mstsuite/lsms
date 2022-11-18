#
# Created by Franco P. Moitzi
#

add_library(mjson STATIC)

target_compile_options(mjson PRIVATE
        $<$<OR:$<C_COMPILER_ID:AppleClang>,$<C_COMPILER_ID:Clang>,$<C_COMPILER_ID:GNU>>:
        -w>)

target_sources(mjson PRIVATE
        ${PROJECT_SOURCE_DIR}/external/mjson/json.c
        ${PROJECT_SOURCE_DIR}/external/mjson/json.h
        ${PROJECT_SOURCE_DIR}/external/mjson/json_helper.c
        ${PROJECT_SOURCE_DIR}/external/mjson/json_helper.h
        ${PROJECT_SOURCE_DIR}/external/mjson/rstring.c
        ${PROJECT_SOURCE_DIR}/external/mjson/rstring.h
        )

target_include_directories(
        mjson PUBLIC
        ${PROJECT_SOURCE_DIR}/external/mjson/
)
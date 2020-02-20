# Sources ---------------------------------------------------------------------
set(SRC
    json.c
    json_helper.c
    rstring.c
)

# Targets ---------------------------------------------------------------------
add_library(mjson STATIC ${SRC})

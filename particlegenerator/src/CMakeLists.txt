#file(GLOB HDRS *.h)
set(HDRS particle2.h
         particles.h
)
set(SRCS particles.cpp
)

include_directories(../../sceneparser/)

add_library(particlegenerator STATIC ${SRCS} ${HDRS} ${MOC_SRCS})
if ( WIN32 )
    target_link_libraries(particlegenerator Qt5::Core Qt5::WinMain sceneparser)
else ()
    target_link_libraries(particlegenerator Qt5::Core sceneparser)
endif ()

INSTALL(TARGETS particlegenerator
        DESTINATION "lib")
INSTALL(FILES   ${HDRS}
        DESTINATION "include/sphase")

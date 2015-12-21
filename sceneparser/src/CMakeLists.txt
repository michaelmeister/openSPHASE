#file(GLOB HDRS *.h)
set(HDRS global.h
         jsonparser.h
         configurationparser.h
         configuration.h
         vector2.h
         linesegment.h
)
set(SRCS jsonparser.cpp
         configurationparser.cpp
         configuration.cpp
)

add_library(sceneparser STATIC ${SRCS} ${HDRS} ${MOC_SRCS})
if ( WIN32 )
    target_link_libraries(sceneparser Qt5::Core Qt5::WinMain)
else ()
    target_link_libraries(sceneparser Qt5::Core)
endif ()

INSTALL(TARGETS sceneparser
        DESTINATION "lib")
INSTALL(FILES   ${HDRS}
        DESTINATION "include/sphase")

#file(GLOB HDRS *.h)
set(HDRS stable.h
         scene2.h
         unigrid.h
         kernels.h
         sph.h
         neighbourdata.h)
set(SRCS main.cpp sph.cpp kernels.cpp unigrid.cpp scene2.cpp)



include_directories(SYSTEM ../../3rdparty/)
add_definitions(-DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_OMP)

set (QT_LIB Qt5::Core)
if (NOT ${NOGUI})
    set(CMAKE_AUTOMOC TRUE)
    include_directories( ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
    QT5_WRAP_UI(UI_HEADERS mainwindow.ui)
    set(HDRS ${HDRS} qcustomplot.h particlerenderer.h)
    set(SRCS ${SRCS} particlerenderer.cpp mainwindow.cpp sphthread.cpp qcustomplot.cpp )
    QT5_ADD_RESOURCES(UI_RESOURCES icons.qrc)
    set (QT_LIB ${QT_LIB} Qt5::Gui Qt5::Xml Qt5::PrintSupport Qt5::OpenGL)
endif (NOT ${NOGUI})

if ( WIN32 )
    set (QT_LIB ${QT_LIB} Qt5::WinMain)
endif ()

add_executable(sph ${SRCS} ${HDRS} ${UI_SRCS} ${MOC_SRCS} ${UI_HEADERS} ${UI_RESOURCES})
target_link_libraries(sph ${QT_LIB} ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY} ${H5Part} ${HDF5} ${QGLVIEWER} GLU particlegenerator sceneparser)

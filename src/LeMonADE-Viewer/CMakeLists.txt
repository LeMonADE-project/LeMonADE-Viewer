cmake_minimum_required (VERSION 2.6)

add_library( Camera Camera.cpp )
target_link_libraries(Camera ${SDL_LIBRARY} ${OPENGL_LIBRARY})

add_library( LineParser LineParser.cpp )
#add_library( LeMonADeOpenGL LeMonADeOpenGL.cpp )

add_executable(${APPLICATION_NAME} LeMonADEViewerMain.cpp)
target_link_libraries(${APPLICATION_NAME} ${OPENGL_LIBRARIES} ${FLTK_LIBRARIES} ${PNG_LIBRARIES} Camera LineParser LeMonADE)
SET_TARGET_PROPERTIES(${APPLICATION_NAME} PROPERTIES VERSION "${APPLICATION_VERSION_MAJOR}.${APPLICATION_VERSION_MINOR}" OUTPUT_NAME ${APPLICATION_NAME} CLEAN_DIRECT_OUTPUT 1)

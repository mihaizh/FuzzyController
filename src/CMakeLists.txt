set(SRC_FILES
	fuzzy.cpp
	fuzzy_mfs.h
	Matrix.h
	Matrix.cpp
	TfRunner.h
	TfRunner.cpp
	FuzzyController.cpp
	FuzzyController.h)
	
find_package(OpenCV REQUIRED)
include_directories(${OPENCV_INCLUDE_DIR})

add_executable(${PROJECT_NAME} ${SRC_FILES})
target_link_libraries(${PROJECT_NAME} ${OpenCV_LIBS})
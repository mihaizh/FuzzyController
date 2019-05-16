set(SRC_FILES
	fuzzy.cpp
	fuzzy_mfs.h
	Matrix.h
	Matrix.cpp
	TfRunner.h
	TfRunner.cpp
	FuzzyController.cpp
	FuzzyController.h)

add_executable(${PROJECT_NAME} ${SRC_FILES})

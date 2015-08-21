
include_directories(
	${CMAKE_SOURCE_DIR}/BaseLib
	${CMAKE_SOURCE_DIR}/FileIO
	${CMAKE_SOURCE_DIR}/GeoLib
	${CMAKE_SOURCE_DIR}/MeshLib
)

## Gocad SGrid to VTK file converter ##
add_executable(GocadSGridReader GocadSGridReaderMain.cpp)
set_target_properties(GocadSGridReader PROPERTIES FOLDER Utilities)
target_link_libraries(GocadSGridReader
	GeoLib
	MeshLib
	FileIO
	${Boost_LIBRARIES}
)



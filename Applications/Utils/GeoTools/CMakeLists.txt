include_directories(
	${CMAKE_SOURCE_DIR}
	${CMAKE_SOURCE_DIR}/BaseLib
	${CMAKE_SOURCE_DIR}/GeoLib
	${CMAKE_SOURCE_DIR}/FileIO
	${CMAKE_SOURCE_DIR}/MathLib
	${CMAKE_SOURCE_DIR}/MeshLib
)

if(QT4_FOUND)
	add_executable(TriangulatePolyline TriangulatePolyline.cpp)
	target_link_libraries(TriangulatePolyline
		BaseLib
		FileIO
		MathLib
		${QT_LIBRARIES}
	)
	set_target_properties(TriangulatePolyline PROPERTIES FOLDER Utilities)
	install(TARGETS
		TriangulatePolyline
		RUNTIME DESTINATION bin
		COMPONENT Utilities
	)
endif() # QT4_FOUND

add_executable(CreateNeumannBCsOnSurface CreateNeumannBCsOnSurface.cpp)
target_link_libraries(CreateNeumannBCsOnSurface
	BaseLib
	FileIO
	MathLib
	MeshLib
	InSituLib
	${CATALYST_LIBRARIES}
)
ADD_CATALYST_DEPENDENCY(CreateNeumannBCsOnSurface)
set_target_properties(CreateNeumannBCsOnSurface PROPERTIES FOLDER Utilities)


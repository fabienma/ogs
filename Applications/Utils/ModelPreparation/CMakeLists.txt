
INCLUDE_DIRECTORIES(
	${CMAKE_SOURCE_DIR}
	${CMAKE_SOURCE_DIR}/BaseLib
	${CMAKE_SOURCE_DIR}/GeoLib
	${CMAKE_SOURCE_DIR}/FileIO
	${CMAKE_SOURCE_DIR}/FileIO/GmshIO
	${CMAKE_SOURCE_DIR}/MathLib
	${CMAKE_SOURCE_DIR}/MeshLib
)

# Create executables
IF(QT4_FOUND)
	ADD_EXECUTABLE( ConvertPolylineToVerticalPolygon
		ConvertPolylineToVerticalPolygon.cpp
		ExtractMeshNodes.cpp
		ExtractMeshNodes.h
	)

	TARGET_LINK_LIBRARIES( ConvertPolylineToVerticalPolygon
		FileIO
		GeoLib
		MeshLib
		BaseLib
		InSituLib
		${CATALYST_LIBRARIES}
		${QT_LIBRARIES}
	)
	ADD_CATALYST_DEPENDENCY(ConvertPolylineToVerticalPolygon)
	set_target_properties(ConvertPolylineToVerticalPolygon PROPERTIES FOLDER Utilities)

ENDIF() # QT4_FOUND



INCLUDE_DIRECTORIES(
	${CMAKE_SOURCE_DIR}
	${CMAKE_SOURCE_DIR}/BaseLib
	${CMAKE_SOURCE_DIR}/GeoLib
	${CMAKE_SOURCE_DIR}/FileIO
	${CMAKE_SOURCE_DIR}/FileIO/GmshIO
	${CMAKE_SOURCE_DIR}/MathLib
	${CMAKE_SOURCE_DIR}/MeshLib
	${CMAKE_SOURCE_DIR}/FemLib
)

# Create executables
IF(QT4_FOUND)

	ADD_EXECUTABLE( removeMeshNodes removeMeshNodes.cpp )
	TARGET_LINK_LIBRARIES( removeMeshNodes
		BaseLib
		FileIO
		MeshLib
		${QT_LIBRARIES}
	)

	ADD_EXECUTABLE( moveMeshNodes moveMeshNodes.cpp )
	TARGET_LINK_LIBRARIES( moveMeshNodes
		BaseLib
		FileIO
		MathLib
		MeshLib
		${QT_LIBRARIES}
	)

	SET_TARGET_PROPERTIES(moveMeshNodes removeMeshNodes
		PROPERTIES FOLDER Utilities)


	ADD_EXECUTABLE( convertPolylineToVerticalPolygon 
		convertPolylineToVerticalPolygon.cpp 
		ExtractMeshNodes.cpp
		ExtractMeshNodes.h
	)
	
	TARGET_LINK_LIBRARIES( convertPolylineToVerticalPolygon
		FileIO
		GeoLib
		MeshLib
		BaseLib
		OgsLib
		${QT_LIBRARIES}
	)

ENDIF() # QT4_FOUND

ADD_EXECUTABLE( DataExplorer5NodeReordering DataExplorer5NodeReordering.cpp )
TARGET_LINK_LIBRARIES( DataExplorer5NodeReordering
	FileIO
	MeshLib
)
SET_TARGET_PROPERTIES(DataExplorer5NodeReordering PROPERTIES FOLDER Utilities)

ADD_EXECUTABLE( MoveMesh MoveMesh.cpp )
TARGET_LINK_LIBRARIES( MoveMesh
	FileIO
	MeshLib
)
SET_TARGET_PROPERTIES(MoveMesh PROPERTIES FOLDER Utilities)


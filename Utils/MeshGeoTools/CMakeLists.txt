INCLUDE_DIRECTORIES(
	${CMAKE_SOURCE_DIR}
	${CMAKE_SOURCE_DIR}/OGS
	${CMAKE_SOURCE_DIR}/BaseLib
	${CMAKE_SOURCE_DIR}/FileIO
	${CMAKE_SOURCE_DIR}/FileIO/Legacy
	${CMAKE_SOURCE_DIR}/GeoLib
	${CMAKE_SOURCE_DIR}/MathLib
	${CMAKE_SOURCE_DIR}/MeshLib
)

ADD_EXECUTABLE( ComputeSurfaceNodeIDsInPolygonalRegion ComputeSurfaceNodeIDsInPolygonalRegion.cpp )
SET_TARGET_PROPERTIES( ComputeSurfaceNodeIDsInPolygonalRegion PROPERTIES FOLDER Utils)
TARGET_LINK_LIBRARIES( ComputeSurfaceNodeIDsInPolygonalRegion
	MeshLib
	FileIO
	GeoLib
	MathLib
	OgsLib
	${Boost_LIBRARIES}
	QtDataView
	VtkVis
	${OGS_VTK_REQUIRED_LIBS}
	${QT_LIBRARIES}
)


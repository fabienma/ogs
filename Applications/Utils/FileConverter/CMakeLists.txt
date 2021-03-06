
include_directories(
	${CMAKE_SOURCE_DIR}/Utils/FileConverter
	${CMAKE_SOURCE_DIR}/BaseLib
	${CMAKE_SOURCE_DIR}/GeoLib
	${CMAKE_SOURCE_DIR}/FileIO
	${CMAKE_SOURCE_DIR}/MathLib
	${CMAKE_SOURCE_DIR}/MeshLib
	${CMAKE_SOURCE_DIR}/Gui/OpenSG
)

if (QT4_FOUND)
	add_executable(ConvertSHPToGLI ConvertSHPToGLI.cpp)
	set_target_properties(ConvertSHPToGLI PROPERTIES FOLDER Utilities)

	target_link_libraries(ConvertSHPToGLI FileIO)

	add_executable(generateBCandGLI generateBCandGLI.cpp )
	target_link_libraries(generateBCandGLI FileIO)

	add_executable(generateBCFromPolyline generateBCFromPolyline.cpp )
	target_link_libraries(generateBCFromPolyline FileIO)

	set_target_properties(generateBCandGLI generateBCFromPolyline
		PROPERTIES FOLDER Utilities)

	add_executable(FEFLOW2OGS FEFLOW2OGS.cpp)
	set_target_properties(FEFLOW2OGS PROPERTIES FOLDER Utilities)
	target_link_libraries(FEFLOW2OGS FileIO)
	ADD_CATALYST_DEPENDENCY(FEFLOW2OGS)
endif () # QT4_FOUND

add_executable(generateMatPropsFromMatID generateMatPropsFromMatID.cpp )
target_link_libraries(generateMatPropsFromMatID FileIO)
ADD_CATALYST_DEPENDENCY(generateMatPropsFromMatID)
set_target_properties(generateMatPropsFromMatID
	PROPERTIES FOLDER Utilities)

add_executable(GMSH2OGS GMSH2OGS.cpp)
set_target_properties(GMSH2OGS PROPERTIES FOLDER Utilities)
target_link_libraries(GMSH2OGS FileIO )
ADD_CATALYST_DEPENDENCY(GMSH2OGS)

add_executable(OGS2VTK OGS2VTK.cpp)
set_target_properties(OGS2VTK PROPERTIES FOLDER Utilities)
target_link_libraries(OGS2VTK FileIO)
ADD_CATALYST_DEPENDENCY(OGS2VTK)

add_executable(VTK2OGS VTK2OGS.cpp)
set_target_properties(VTK2OGS PROPERTIES FOLDER Utilities)
target_link_libraries(VTK2OGS FileIO)
ADD_CATALYST_DEPENDENCY(VTK2OGS)

add_executable(VTK2TIN VTK2TIN.cpp)
set_target_properties(VTK2TIN PROPERTIES FOLDER Utilities)
target_link_libraries(VTK2TIN FileIO)
ADD_CATALYST_DEPENDENCY(VTK2TIN)

add_executable(TIN2VTK TIN2VTK.cpp)
set_target_properties(TIN2VTK PROPERTIES FOLDER Utilities)
target_link_libraries(TIN2VTK FileIO)
ADD_CATALYST_DEPENDENCY(TIN2VTK)

####################
### Installation ###
####################
install(TARGETS generateMatPropsFromMatID GMSH2OGS OGS2VTK VTK2OGS VTK2TIN
	RUNTIME DESTINATION bin COMPONENT ogs_converter)

if(QT4_FOUND)
	install(TARGETS ConvertSHPToGLI generateBCandGLI generateBCFromPolyline generateBCandGLI
		FEFLOW2OGS RUNTIME DESTINATION bin COMPONENT ogs_converter)
endif()

cpack_add_component(ogs_converter
	DISPLAY_NAME "File converter"
	DESCRIPTION "File converter."
	GROUP Utilities
)

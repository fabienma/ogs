include_directories(
	${CMAKE_SOURCE_DIR}
	${CMAKE_SOURCE_DIR}/BaseLib
	${CMAKE_SOURCE_DIR}/GeoLib
	${CMAKE_SOURCE_DIR}/FileIO
	${CMAKE_SOURCE_DIR}/MathLib
)

if(QT4_FOUND)
	add_executable(CreateTransientNeumannConditions
	CreateTransientNeumannConditions.cpp)
	target_link_libraries(CreateTransientNeumannConditions
		BaseLib
		FileIO
		MathLib
		InSituLib
		${CATALYST_LIBRARIES}
	)
	set_target_properties(CreateTransientNeumannConditions PROPERTIES FOLDER Utilities)
	install(TARGETS
		CreateTransientNeumannConditions
		RUNTIME DESTINATION bin
		COMPONENT Utilities
	)
endif() # QT4_FOUND

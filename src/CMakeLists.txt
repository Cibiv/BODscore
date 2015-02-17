cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -std=c++0x")

include_directories (sambamio)

add_executable(vshape 	BasisClass.cpp
                        Coverage.cpp
						Ident_hist.cpp
						ParseSNP.cpp
						nvshape.cpp     #VShape.cpp
						sambamio/Alignment.cpp
						sambamio/BamParser.cpp
						sambamio/FastaParser.cpp
						sambamio/Parser.cpp
						sambamio/SamParser.cpp
 					)
 					
TARGET_LINK_LIBRARIES(vshape BamTools-static)
TARGET_LINK_LIBRARIES(vshape zlibstatic)
target_link_libraries(vshape sqlite3pp)

install (TARGETS vshape
         RUNTIME DESTINATION /usr/bin/
         LIBRARY DESTINATION /usr/lib/
         ARCHIVE DESTINATION /usr/lib/)



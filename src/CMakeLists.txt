cmake_minimum_required(VERSION 2.8)

include_directories (sambamio)

add_executable(vshape 	BasisClass.cpp
						Ident_hist.cpp
						ParseSNP.cpp
						VShape.cpp
						sambamio/Alignment.cpp
						sambamio/BamParser.cpp
						sambamio/FastaParser.cpp
						sambamio/Parser.cpp
						sambamio/SamParser.cpp
 					)
 					
TARGET_LINK_LIBRARIES(vshape BamTools-static)
TARGET_LINK_LIBRARIES(vshape zlibstatic)
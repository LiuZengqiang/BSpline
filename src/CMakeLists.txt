message('-----BEGIN B_SPLINE_LIB-----')

add_library(BSpline_lib SHARED BSpline.cpp)

target_link_libraries(BSpline_lib PUBLIC ${THIRD_PARTY_LIBS})

message('-----END B_SPLINE_LIB-----')
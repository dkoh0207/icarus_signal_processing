add_subdirectory(Filters)
add_subdirectory(Detection)

add_definitions(-DEIGEN_FFTW_DEFAULT)

art_make( 
          LIB_LIBRARIES 
                          icarus_signal_processing_Detection
                          icarus_signal_processing_Filters
                          ${ICARUS_FFTW_LIBRARIES}
                          ${ROOT_GEOM}
                          ${ROOT_XMLIO}
                          ${ROOT_GDML}
                          ${ROOT_FFTW}
                          ${ROOT_BASIC_LIB_LIST}
       )

install_headers()
install_fhicl()
install_source()

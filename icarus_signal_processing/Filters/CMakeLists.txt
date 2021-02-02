add_definitions(-DEIGEN_FFTW_DEFAULT)

art_make( 
          LIB_LIBRARIES 
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

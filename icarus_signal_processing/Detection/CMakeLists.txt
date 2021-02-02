
art_make( 
          LIB_LIBRARIES 
                          icarus_signal_processing_Filters
                          ${ROOT_GEOM}
                          ${ROOT_XMLIO}
                          ${ROOT_GDML}
                          ${ROOT_FFTW}
                          ${ROOT_BASIC_LIB_LIST}
       )

install_headers()
install_fhicl()
install_source()

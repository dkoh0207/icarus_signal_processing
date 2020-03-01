//
// cint script to generate libraries
// Declaire namespace & classes you defined
// #pragma statement: order matters! Google it ;)
//

#ifdef __CLING__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class sigproc_tools::Morph1D+;
#pragma link C++ class sigproc_tools::Morph2D+;
#pragma link C++ class sigproc_tools::MiscUtils+;
#pragma link C++ class sigproc_tools::Denoising+;
#pragma link C++ class sigproc_tools::AdaptiveWiener+;
#pragma link C++ class Deconvolution::sigproc_tools+;
//ADD_NEW_CLASS ... do not change this line
#endif


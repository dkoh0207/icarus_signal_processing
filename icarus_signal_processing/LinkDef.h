//
// cint script to generate libraries
// Declaire namespace & classes you defined
// #pragma statement: order matters! Google it ;)
//

#ifdef __CLING__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class icarus_signal_processing::Denoising+;
#pragma link C++ class icarus_signal_processing::FindROI2D+;
#pragma link C++ class ICARUSSigProcDefs+;
#pragma link C++ class icarus_signal_processing::WaveformTools+;
//ADD_NEW_CLASS ... do not change this line
#endif


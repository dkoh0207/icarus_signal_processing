//
// cint script to generate libraries
// Declaire namespace & classes you defined
// #pragma statement: order matters! Google it ;)
//

#ifdef __CLING__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class icarus_signal_processing::Filters::AdaptiveWiener+;
#pragma link C++ class icarus_signal_processing::Filters::BilateralFilters+;
#pragma link C++ class icarus_signal_processing::Filters::FFTFilterFunctions+;
#pragma link C++ class icarus_signal_processing::Filters::ICARUSFFT+;
#pragma link C++ class icarus_signal_processing::Filters::ImageFilters+;
#pragma link C++ class icarus_signal_processing::Filters::MiscUtils+;
//ADD_NEW_CLASS ... do not change this line
#endif


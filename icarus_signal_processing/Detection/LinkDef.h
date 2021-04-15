//
// cint script to generate libraries
// Declaire namespace & classes you defined
// #pragma statement: order matters! Google it ;)
//

#ifdef __CLING__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class icarus_signal_processing::Detection::EdgeDetection+;
#pragma link C++ class icarus_signal_processing::Detection::LineDetection+;
#pragma link C++ class icarus_signal_processing::Detection::MorphologicalFunctions1D+;
#pragma link C++ class icarus_signal_processing::Detection::MorphologicalFunctions2D+;
#pragma link C++ class icarus_signal_processing::Detection::Thresholding+;
#pragma link C++ class icarus_signal_processing::Detection::DisjointSetForest+;
//ADD_NEW_CLASS ... do not change this line
#endif


#include "TLorentzVector.h"
#include "vector"
#pragma clang diagnostic ignored "-Wshadow-ivar"
#ifdef __CINT__ 
#pragma link C++ nestedclasses;
#pragma link C++ nestedtypedefs;
#pragma link C++ class vector<TLorentzVector>+;
// #pragma link C++ class vector<TLorentzVector>::*;
#pragma link C++ class vector<vector<vector<Float_t>>>+;
// #pragma link C++ class vector<vector<vector<Float_t>>>::*;
#ifdef G__VECTOR_HAS_CLASS_ITERATOR
#pragma link C++ operators vector<TLorentzVector>::iterator;
#pragma link C++ operators vector<TLorentzVector>::const_iterator;
#pragma link C++ operators vector<TLorentzVector>::reverse_iterator;
#pragma link C++ operators vector<vector<vector<Float_t>>>::iterator;
#pragma link C++ operators vector<vector<vector<Float_t>>>::const_iterator;
#pragma link C++ operators vector<vector<vector<Float_t>>>::reverse_iterator;
#endif
#endif

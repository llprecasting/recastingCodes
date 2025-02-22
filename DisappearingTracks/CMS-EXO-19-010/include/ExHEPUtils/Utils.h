// -*- C++ -*-
//
// This file is part of HEPUtils -- https://bitbucket.org/andybuckley/heputils
// Copyright (C) 2013-2018 Andy Buckley <andy.buckley@cern.ch>
//
// Embedding of HEPUtils code in other projects is permitted provided this
// notice is retained and the HEPUtils namespace and include path are changed.
//

// Duly: HEPUtils namespace changed to HEP

#pragma once

#if defined(__cplusplus) && __cplusplus < 201103L
// #define XSTR(x) STR(x)
// #define STR(x) #x
// #pragma message STR(__cplusplus)
#error "This library needs at least a C++11 compliant compiler"
#endif

/// @file Utility functions
/// @author Andy Buckley <andy.buckley@cern.ch>

namespace HEP {


  /// @name Container utils
  //@{

  /// Return true if f(x) is true for any x in container c, otherwise false.
  template <typename CONTAINER, typename FN>
  inline bool any(const CONTAINER& c, const FN& f) {
    for (const typename CONTAINER::value_type& val : c) {
      if (f(val)) return true;
    }
    return false;
  }

  /// Return true if f(x) is true for all x in container c, otherwise false.
  template <typename CONTAINER, typename FN>
  inline bool all(const CONTAINER& c, const FN& f) {
    for (const typename CONTAINER::value_type& val : c) {
      if (!f(val)) return false;
    }
    return true;
  }

  /// Return true if x is found in container c, otherwise false.
  template <typename CONTAINER, typename T>
  inline bool contains(const CONTAINER& c, const T& x) {
    for (const typename CONTAINER::value_type& val : c) {
      if (val == x) return true;
    }
    return false;
  }


  //@}


}

#ifndef AMREX_EB2_IF_SPLINE_H_
#define AMREX_EB2_IF_SPLINE_H_

#include <AMReX_Array.H>
#include "AMReX_distFcnElement.H"
#include <algorithm>
#include "AMReX_Vector.H"

// For all implicit functions, >0: body; =0: boundary; <0: fluid

namespace amrex { namespace EB2 {

class SplineIF
{
public:
    // inside: is the fluid inside the spline element?

  SplineIF (std::string splinefile):m_spline_file(splinefile)
  {
    // read_spline_profile(m_spline_file, m_splpts);
    // addSplineElement(splpts);
  }

  ~SplineIF () {
    for (int i=0; i<geomElements.size(); ++i) {
      delete geomElements[i];
    }
    geomElements.clear();
  }
  
  // SplineIF (const SplineIF& rhs) noexcept = default;
  // SplineIF (SplineIF&& rhs) noexcept = default;
  // SplineIF& operator= (const SplineIF& rhs) = delete;
  // SplineIF& operator= (SplineIF&& rhs) = delete;
  
  Real operator() (const RealArray& p) const
  {
    amrex::RealVect cp;
    distFcnElement2d * closesetGeomElement;
    amrex::Real dist;
    dist = 1.0e29;
    for (auto * geom : geomElements ) {
      amrex::Real d = geom->cpdist(p, cp);
      if (d < dist) {
        dist = d;
        closesetGeomElement = geom;
      }
    }
    amrex::Real side = closesetGeomElement->cpside(p,cp);
    return dist*side;
  }

  void addGeomElement(distFcnElement2d * newElement) {
    geomElements.push_back(newElement);
  }
  void addLineElement(std::vector<amrex::RealVect> pts){
    LineDistFcnElement2d * theLine = new LineDistFcnElement2d();
    theLine->set_control_points(pts);
    geomElements.push_back(theLine);
  }
  void addSplineElement(std::vector<amrex::RealVect> pts){
    SplineDistFcnElement2d * theSpline = new SplineDistFcnElement2d();
    theSpline->set_control_points(pts);
    theSpline->calc_D();
    geomElements.push_back(theSpline);
  }

protected:

  std::string m_spline_file;
  RealArray m_splpts;
  // The geometry elements used to compute distance function
  amrex::Vector<distFcnElement2d*> geomElements;
};

}}

#endif

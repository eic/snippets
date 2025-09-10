/***************************************************************************
 *
 * $Id: StPhysicalHelix.hh,v 1.4 2005/07/06 18:49:56 fisyak Exp $
 *
 * Author: Brian Lasiuk, Sep 1997
 ***************************************************************************
 *
 * Description: 
 * Parametrization of a physical helix. See the SCL user guide for more.
 * 
 ***************************************************************************
 *
 * $Log: StPhysicalHelix.hh,v $
 * Revision 1.4  2005/07/06 18:49:56  fisyak
 * Replace StHelixD, StLorentzVectorD,StLorentzVectorF,StMatrixD,StMatrixF,StPhysicalHelixD,StThreeVectorD,StThreeVectorF by templated version
 *
 * Revision 1.3  2002/06/21 17:49:26  genevb
 * Some minor speed improvements
 *
 * Revision 1.2  2002/02/20 00:56:23  ullrich
 * Added methods to calculate signed DCA.
 *
 * Revision 1.1  1999/01/30 03:59:04  fisyak
 * Root Version of StarClassLibrary
 *
 * Revision 1.1  1999/01/23 00:27:59  ullrich
 * Initial Revision
 *
 **************************************************************************/
#ifndef ST_PHYSICAL_HELIX_HH
#define ST_PHYSICAL_HELIX_HH

#include "TVector3.h"
#include "StHelix.h"

class StPhysicalHelix : public StHelix {
public:
    // Requires: momentum, origin, signed Magnetic Field
    //           and Charge of particle (+/- 1)
    StPhysicalHelix(const TVector3&,
		    const TVector3&,
		    double, double);
    
    // curvature, dip angle, phase, origin, h
    StPhysicalHelix(double, double, double,
		    const TVector3&, int h=-1);
    StPhysicalHelix();
    
    ~StPhysicalHelix();

    // Requires:  signed Magnetic Field
    TVector3 momentum(double) const;     // returns the momentum at origin
    TVector3 momentumAt(double, double) const; // returns momemtum at S
    int                   charge(double)   const;     // returns charge of particle
    // 2d DCA to x,y point signed relative to curvature
    double curvatureSignedDistance(double x, double y) ;
    // 2d DCA to x,y point signed relative to rotation 
    double geometricSignedDistance(double x, double y) ;
    // 3d DCA to 3d point signed relative to curvature
    double curvatureSignedDistance(const TVector3&) ;
    // 3d DCA to 3d point signed relative to rotation
    double geometricSignedDistance(const TVector3&) ;
    
#ifdef __ROOT__
  ClassDef(StPhysicalHelix,1)
#endif
};

#endif

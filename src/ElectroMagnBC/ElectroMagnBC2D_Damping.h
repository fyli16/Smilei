
#ifndef ELECTROMAGNBC2D_DAMPING_H
#define ELECTROMAGNBC2D_DAMPING_H

#include "ElectroMagnBC2D.h"

class Params;
class ElectroMagn;

//! Field-damping "mask" boundary condition (2D Cartesian).
//!
//! Implements an absorbing sponge layer of finite thickness placed *inside* the
//! simulation box, at x-min (i_boundary=0) or x-max (i_boundary=1). Inside the layer
//! the transverse electromagnetic wave fields (Ey, Ez, By, Bz) are multiplied every
//! time step by a smooth coefficient that ramps from 1 (no damping) at the inner edge
//! of the mask to (1 - coefficient) at the box wall. This progressively absorbs the
//! Alfven wave before it reaches the box boundary.
//!
//! The parallel guide field Bx (uniform B0) and the parallel electric field Ex are
//! left untouched, so the mean field and the electrostatic plasma response are preserved.
//! The mask only depends on x, so the coefficient is uniform along y.
class ElectroMagnBC2D_Damping : public ElectroMagnBC2D
{
public:
    ElectroMagnBC2D_Damping( Params &params, Patch *patch, unsigned int i_boundary );
    ~ElectroMagnBC2D_Damping() {};

    void apply( ElectroMagn *EMfields, double time_dual, Patch *patch ) override;

private:

    //! Damping coefficient at a given physical position x (returns 1 outside the mask)
    double dampingCoefficient( double x ) const;

    //! Thickness of the mask layer (in code length units)
    double thickness_;

    //! Maximum damping strength in [0,1] (1 => field fully removed at the wall)
    double coefficient_;

    //! Upper bound of the box along x (code length units)
    double x_max_;

    //! Cell length along x
    double dx_;

    //! Number of ghost cells along x
    unsigned int oversize_x_;
};

#endif


#include "ElectroMagnBC2D_Damping.h"

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>

#include "Params.h"
#include "Patch.h"
#include "ElectroMagn2D.h"
#include "Field2D.h"
#include "Tools.h"

using namespace std;

// ---------------------------------------------------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------------------------------------------------
ElectroMagnBC2D_Damping::ElectroMagnBC2D_Damping( Params &params, Patch *patch, unsigned int i_boundary )
    : ElectroMagnBC2D( params, patch, i_boundary )
{
    thickness_   = params.EM_damping_thickness;
    coefficient_ = params.EM_damping_coefficient;
    x_max_       = params.grid_length[0];
    dx_          = d[0];
    oversize_x_  = params.oversize[0];

    if( thickness_ <= 0. ) {
        ERROR( "EM boundary condition 'damping' requires a strictly positive `EM_damping_thickness`" );
    }
    if( coefficient_ < 0. || coefficient_ > 1. ) {
        ERROR( "`EM_damping_coefficient` must lie in [0,1] (got " << coefficient_ << ")" );
    }
    if( 2.*thickness_ >= x_max_ ) {
        ERROR( "`EM_damping_thickness` (" << thickness_ << ") is too large: the two masks would overlap "
               << "(box length along x = " << x_max_ << ")" );
    }
}

// ---------------------------------------------------------------------------------------------------------------------
// Smooth parabolic damping profile, = 1 outside the mask, ramping to (1-coefficient_) at the wall
// ---------------------------------------------------------------------------------------------------------------------
double ElectroMagnBC2D_Damping::dampingCoefficient( double x ) const
{
    double s; // normalized depth into the mask: 0 at the inner edge, 1 at the box wall
    if( i_boundary_ == 0 ) {
        // Left mask: [0, thickness_]
        if( x >= thickness_ ) {
            return 1.;
        }
        s = ( thickness_ - x ) / thickness_;
    } else {
        // Right mask: [x_max_ - thickness_, x_max_]
        if( x <= x_max_ - thickness_ ) {
            return 1.;
        }
        s = ( x - ( x_max_ - thickness_ ) ) / thickness_;
    }
    if( s < 0. ) {
        s = 0.;
    } else if( s > 1. ) {
        s = 1.;
    }
    return 1. - coefficient_ * s * s;
}

// ---------------------------------------------------------------------------------------------------------------------
// Apply the field-damping mask
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnBC2D_Damping::apply( ElectroMagn *EMfields, double, Patch *patch )
{
    // Physical position (along x) of a node, given its local x-index and grid staggering.
    // shift = 0.0 for a node primal in x, 0.5 for a node dual in x.
    const double x_patch_min = patch->getDomainLocalMin( 0 );
    auto nodePosition = [&]( unsigned int i, double shift ) -> double {
        return x_patch_min + ( ( double )i - ( double )oversize_x_ + shift ) * dx_;
    };

    // Transverse wave fields only. The guide field Bx (uniform B0) and the parallel
    // field Ex are intentionally left untouched.
    Field2D *Ey2D = static_cast<Field2D *>( EMfields->Ey_ ); // primal in x, dual   in y
    Field2D *Ez2D = static_cast<Field2D *>( EMfields->Ez_ ); // primal in x, primal in y
    Field2D *By2D = static_cast<Field2D *>( EMfields->By_ ); // dual   in x, primal in y
    Field2D *Bz2D = static_cast<Field2D *>( EMfields->Bz_ ); // dual   in x, dual   in y

    // Primal-in-x fields: Ey (x:primal, y:dual), Ez (x:primal, y:primal)
    for( unsigned int i=0 ; i<n_p[0] ; i++ ) {
        const double c = dampingCoefficient( nodePosition( i, 0. ) );
        if( c == 1. ) {
            continue;
        }
        for( unsigned int j=0 ; j<n_d[1] ; j++ ) {
            ( *Ey2D )( i, j ) *= c;
        }
        for( unsigned int j=0 ; j<n_p[1] ; j++ ) {
            ( *Ez2D )( i, j ) *= c;
        }
    }

    // Dual-in-x fields: By (x:dual, y:primal), Bz (x:dual, y:dual)
    for( unsigned int i=0 ; i<n_d[0] ; i++ ) {
        const double c = dampingCoefficient( nodePosition( i, 0.5 ) );
        if( c == 1. ) {
            continue;
        }
        for( unsigned int j=0 ; j<n_p[1] ; j++ ) {
            ( *By2D )( i, j ) *= c;
        }
        for( unsigned int j=0 ; j<n_d[1] ; j++ ) {
            ( *Bz2D )( i, j ) *= c;
        }
    }
}

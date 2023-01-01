// ----------------------------------------------------------------------------
//      ______  __  __  ______
//     / ____/ / / / / / ____/
//    / / __  / /_/ / / /___
//   / / / / / __  / /___  /
//  / /_/ / / / / / ____/ /    GeneralisedHyperbolicStretch
//  \____/ /_/ /_/ /_____/     https://www.ghsastro.co.uk
// ----------------------------------------------------------------------------
// Standard GHS Process Module
// ----------------------------------------------------------------------------
// GHSParameters.cpp
// ----------------------------------------------------------------------------
//
// Copyright (C) 2022,2023 Mike Cranfield
//
// This product is based on software from the PixInsight project, developed
// by Pleiades Astrophoto and its contributors (https://pixinsight.com/).
// ----------------------------------------------------------------------------
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation, version 3 of the License.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License along with
// this program.  If not, see <http://www.gnu.org/licenses/>.
// ----------------------------------------------------------------------------

#include "GHSParameters.h"

namespace pcl
{

// ----------------------------------------------------------------------------

GHSST*                TheGHSSTParameter = nullptr;
GHSSC*                TheGHSSCParameter = nullptr;
GHSInv*               TheGHSInvParameter = nullptr;
GHS*                  TheGHSParameter = nullptr;
GHSb*                 TheGHSbParameter = nullptr;
GHSSP*                TheGHSSPParameter = nullptr;
GHSHP*                TheGHSHPParameter = nullptr;
GHSLP*                TheGHSLPParameter = nullptr;
GHSBP*                TheGHSBPParameter = nullptr;
GHSWP*                TheGHSWPParameter = nullptr;
GHSCB*                TheGHSCBParameter = nullptr;
GHSCT*                TheGHSCTParameter = nullptr;
GHSRGBWS*             TheGHSRGBWSParameter = nullptr;

// ----------------------------------------------------------------------------

GHSST::GHSST( MetaProcess* P ) : MetaEnumeration( P )
{
   TheGHSSTParameter = this;
}

IsoString GHSST::Id() const
{
   return "stretchType";
}

size_type GHSST::NumberOfElements() const
{
   return ST_NumberOfItems;
}

IsoString GHSST::ElementId( size_type i ) const
{
    switch ( i )
    {
        default:
        case ST_GeneralisedHyperbolic:  return "ST_GeneralisedHyperbolic";
        case ST_MidtonesTransfer:       return "ST_MidtonesTransfer";
        case ST_Arcsinh:                return "ST_Arcsinh";
        case ST_Linear:                 return "ST_Linear";
    }
}

int GHSST::ElementValue( size_type i ) const
{
   return int( i );
}

size_type GHSST::DefaultValueIndex() const
{
   return Default;
}

// ----------------------------------------------------------------------------

GHSSC::GHSSC( MetaProcess* P ) : MetaEnumeration( P )
{
    TheGHSSCParameter = this;
}

IsoString GHSSC::Id() const
{
   return "stretchChannel";
}

size_type GHSSC::NumberOfElements() const
{
    return SC_NumberOfItems;
}

IsoString GHSSC::ElementId( size_type i ) const
{
    switch ( i )
    {
        case SC_Red:        return "SC_Red";
        case SC_Green:      return "SC_Green";
        case SC_Blue:       return "SC_Blue";
        case SC_RGB:        return "SC_RGB";
        case SC_Lightness:  return "SC_Lightness";
        case SC_Saturation: return "SC_Saturation";
        case SC_Colour:     return "SC_Colour";
        default:            return "SC_RGB";
        
    }
}

int GHSSC::ElementValue( size_type i ) const
{
   return int( i );
}

size_type GHSSC::DefaultValueIndex() const
{
   return Default;
}

// ----------------------------------------------------------------------------

GHSInv::GHSInv( MetaProcess* P ) : MetaBoolean( P )
{
   TheGHSInvParameter = this;
}

IsoString GHSInv::Id() const
{
   return "inverse";
}

bool GHSInv::DefaultValue() const
{
   return false;
}

// ----------------------------------------------------------------------------

GHS::GHS( MetaProcess* P ) : MetaDouble( P )
{
   TheGHSParameter = this;
}

IsoString GHS::Id() const
{
   return "stretchFactor";
}

int GHS::Precision() const
{
   return 3;
}

double GHS::DefaultValue() const
{
   return 0.0;
}

double GHS::MinimumValue() const
{
   return 0.0;
}

double GHS::MaximumValue() const
{
   return 20.0;
}

// ----------------------------------------------------------------------------

GHSb::GHSb( MetaProcess* P ) : MetaDouble( P )
{
   TheGHSbParameter = this;
}

IsoString GHSb::Id() const
{
   return "localIntensity";
}

int GHSb::Precision() const
{
   return 3;
}

double GHSb::DefaultValue() const
{
   return 0.0;
}

double GHSb::MinimumValue() const
{
   return -5.0;
}

double GHSb::MaximumValue() const
{
   return 15.0;
}

// ----------------------------------------------------------------------------

GHSSP::GHSSP( MetaProcess* P ) : MetaDouble( P )
{
   TheGHSSPParameter = this;
}

IsoString GHSSP::Id() const
{
   return "symmetryPoint";
}

int GHSSP::Precision() const
{
   return 6;
}

double GHSSP::DefaultValue() const
{
   return 0.0;
}

double GHSSP::MinimumValue() const
{
   return 0.0;
}

double GHSSP::MaximumValue() const
{
   return 1.0;
}

// ----------------------------------------------------------------------------

GHSLP::GHSLP( MetaProcess* P ) : MetaDouble( P )
{
   TheGHSLPParameter = this;
}

IsoString GHSLP::Id() const
{
   return "shadowProtection";
}

int GHSLP::Precision() const
{
   return 6;
}

double GHSLP::DefaultValue() const
{
   return 0.0;
}

double GHSLP::MinimumValue() const
{
   return 0.0;
}

double GHSLP::MaximumValue() const
{
   return 1.0;
}

// ----------------------------------------------------------------------------

GHSHP::GHSHP( MetaProcess* P ) : MetaDouble( P )
{
   TheGHSHPParameter = this;
}

IsoString GHSHP::Id() const
{
   return "highlightProtection";
}

int GHSHP::Precision() const
{
   return 6;
}

double GHSHP::DefaultValue() const
{
   return 1.0;
}

double GHSHP::MinimumValue() const
{
   return 0.0;
}

double GHSHP::MaximumValue() const
{
   return 1.0;
}

// ----------------------------------------------------------------------------

GHSBP::GHSBP( MetaProcess* P ) : MetaDouble( P )
{
   TheGHSBPParameter = this;
}

IsoString GHSBP::Id() const
{
   return "blackPoint";
}

int GHSBP::Precision() const
{
   return 6;
}

double GHSBP::DefaultValue() const
{
   return 0.0;
}

double GHSBP::MinimumValue() const
{
    return -1.0;
}

double GHSBP::MaximumValue() const
{
   return 1.0;
}

// ----------------------------------------------------------------------------

GHSWP::GHSWP( MetaProcess* P ) : MetaDouble( P )
{
   TheGHSWPParameter = this;
}

IsoString GHSWP::Id() const
{
   return "whitePoint";
}

int GHSWP::Precision() const
{
   return 6;
}

double GHSWP::DefaultValue() const
{
   return 1.0;
}

double GHSWP::MinimumValue() const
{
    return 0.0;
}

double GHSWP::MaximumValue() const
{
   return 2.0;
}


// ----------------------------------------------------------------------------

GHSCB::GHSCB( MetaProcess* P ) : MetaDouble( P )
{
   TheGHSCBParameter = this;
}

IsoString GHSCB::Id() const
{
   return "colourBlend";
}

int GHSCB::Precision() const
{
   return 3;
}

double GHSCB::DefaultValue() const
{
   return 1.0;
}

double GHSCB::MinimumValue() const
{
    return 0.0;
}

double GHSCB::MaximumValue() const
{
   return 1.0;
}

// ----------------------------------------------------------------------------

GHSCT::GHSCT( MetaProcess* P ) : MetaEnumeration( P )
{
    TheGHSCTParameter = this;
}

IsoString GHSCT::Id() const
{
   return "clipType";
}

size_type GHSCT::NumberOfElements() const
{
   return CT_NumberOfItems;
}

IsoString GHSCT::ElementId( size_type i ) const
{
    switch ( i )
    {
        case CT_Clip:           return "CT_Clip";
        case CT_Rescale:        return "CT_Rescale";
        case CT_RGBBlend:       return "CT_RGBBlend";
        case CT_RescaleGlobal:  return "CT_RescaleGlobal";
        default:                return "CT_Clip";
    }
}

int GHSCT::ElementValue( size_type i ) const
{
   return int( i );
}

size_type GHSCT::DefaultValueIndex() const
{
   return Default;
}

// ----------------------------------------------------------------------------

GHSRGBWS::GHSRGBWS( MetaProcess* P ) : MetaBoolean( P )
{
   TheGHSRGBWSParameter = this;
}

IsoString GHSRGBWS::Id() const
{
   return "useRGBWorkingSpace";
}

bool GHSRGBWS::DefaultValue() const
{
   return false;
}

// ----------------------------------------------------------------------------


} // pcl

// ----------------------------------------------------------------------------
// EOF GHSParameters.cpp

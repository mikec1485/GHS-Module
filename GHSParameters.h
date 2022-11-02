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
// GHSParameters.h
// ----------------------------------------------------------------------------
//
// Copyright (C) 2022, Mike Cranfield
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

#ifndef __GHSParameters_h
#define __GHSParameters_h

#include <pcl/MetaParameter.h>

namespace pcl
{

// ----------------------------------------------------------------------------

PCL_BEGIN_LOCAL

// ----------------------------------------------------------------------------

class GHSST : public MetaEnumeration
{
public:

   enum { ST_GeneralisedHyperbolic,
          ST_MidtonesTransfer,
          ST_Arcsinh,
          ST_Linear,
          ST_NumberOfItems,
          Default = ST_GeneralisedHyperbolic };

    GHSST( MetaProcess* );

   IsoString Id() const override;
   size_type NumberOfElements() const override;
   IsoString ElementId( size_type ) const override;
   int ElementValue( size_type ) const override;
   size_type DefaultValueIndex() const override;
};

extern GHSST* TheGHSSTParameter;

// ----------------------------------------------------------------------------

class GHSSC : public MetaEnumeration
{
public:

   enum { SC_Red,
          SC_Green,
          SC_Blue,
          SC_RGB,
          SC_Lightness,
          SC_Saturation,
          SC_Colour,
          SC_NumberOfItems,
          Default = SC_RGB };

    GHSSC( MetaProcess* );

   IsoString Id() const override;
   size_type NumberOfElements() const override;
   IsoString ElementId( size_type ) const override;
   int ElementValue( size_type ) const override;
   size_type DefaultValueIndex() const override;
};

extern GHSSC* TheGHSSCParameter;

// ----------------------------------------------------------------------------

class GHSInv : public MetaBoolean
{
public:

   GHSInv( MetaProcess* );

   IsoString Id() const override;
   bool DefaultValue() const override;
};

extern GHSInv* TheGHSInvParameter;

// ----------------------------------------------------------------------------

class GHS : public MetaDouble
{
public:

   GHS( MetaProcess* );

   IsoString Id() const override;
   int Precision() const override;
   double DefaultValue() const override;
   double MinimumValue() const override;
   double MaximumValue() const override;
};

extern GHS* TheGHSParameter;

// ----------------------------------------------------------------------------

class GHSb : public MetaDouble
{
public:

   GHSb( MetaProcess* );

   IsoString Id() const override;
   int Precision() const override;
   double DefaultValue() const override;
   double MinimumValue() const override;
   double MaximumValue() const override;
};

extern GHSb* TheGHSbParameter;

// ----------------------------------------------------------------------------

class GHSSP : public MetaDouble
{
public:

   GHSSP( MetaProcess* );

   IsoString Id() const override;
   int Precision() const override;
   double DefaultValue() const override;
   double MinimumValue() const override;
   double MaximumValue() const override;
};

extern GHSSP* TheGHSSPParameter;

// ----------------------------------------------------------------------------

class GHSHP : public MetaDouble
{
public:

   GHSHP( MetaProcess* );

   IsoString Id() const override;
   int Precision() const override;
   double DefaultValue() const override;
   double MinimumValue() const override;
   double MaximumValue() const override;
};

extern GHSHP* TheGHSHPParameter;

// ----------------------------------------------------------------------------

class GHSLP : public MetaDouble
{
public:

   GHSLP( MetaProcess* );

   IsoString Id() const override;
   int Precision() const override;
   double DefaultValue() const override;
   double MinimumValue() const override;
   double MaximumValue() const override;
};

extern GHSLP* TheGHSLPParameter;

// ----------------------------------------------------------------------------

class GHSBP : public MetaDouble
{
public:

   GHSBP(MetaProcess*);

   IsoString Id() const override;
   int Precision() const override;
   double DefaultValue() const override;
   double MinimumValue() const override;
   double MaximumValue() const override;
};

extern GHSBP* TheGHSBPParameter;

// ----------------------------------------------------------------------------

class GHSWP : public MetaDouble
{
public:

    GHSWP(MetaProcess*);

   IsoString Id() const override;
   int Precision() const override;
   double DefaultValue() const override;
   double MinimumValue() const override;
   double MaximumValue() const override;
};

extern GHSWP* TheGHSWPParameter;

// ----------------------------------------------------------------------------

class GHSCB : public MetaDouble
{
public:

    GHSCB(MetaProcess*);

   IsoString Id() const override;
   int Precision() const override;
   double DefaultValue() const override;
   double MinimumValue() const override;
   double MaximumValue() const override;
};

extern GHSCB* TheGHSCBParameter;

// ----------------------------------------------------------------------------

class GHSCT : public MetaEnumeration
{
public:

   enum { CT_Clip,
          CT_Rescale,
          CT_RGBBlend,
          CT_RescaleGlobal,
          CT_NumberOfItems,
          Default = CT_Clip };

    GHSCT( MetaProcess* );

   IsoString Id() const override;
   size_type NumberOfElements() const override;
   IsoString ElementId( size_type ) const override;
   int ElementValue( size_type ) const override;
   size_type DefaultValueIndex() const override;
};

extern GHSCT* TheGHSCTParameter;

// ----------------------------------------------------------------------------

class GHSRGBWS : public MetaBoolean
{
public:

   GHSRGBWS( MetaProcess* );

   IsoString Id() const override;
   bool DefaultValue() const override;
};

extern GHSRGBWS* TheGHSRGBWSParameter;

// ----------------------------------------------------------------------------

PCL_END_LOCAL

// ----------------------------------------------------------------------------

} // pcl

#endif   // __GHSParameters_h

// ----------------------------------------------------------------------------
// EOF GHSParameters.h

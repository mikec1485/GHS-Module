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
// GHSProcess.h
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

#ifndef __GHSProcess_h
#define __GHSProcess_h

#include <pcl/MetaProcess.h>

namespace pcl
{

// ----------------------------------------------------------------------------

class GHSProcess : public MetaProcess
{
public:

   GHSProcess();

   IsoString Id() const override;
   IsoString Categories() const override;
   uint32 Version() const override;
   String Description() const override;
   String IconImageSVGFile() const override;
   ProcessInterface* DefaultInterface() const override;
   ProcessImplementation* Create() const override;
   ProcessImplementation* Clone( const ProcessImplementation& ) const override;
   bool CanProcessCommandLines() const override;
   int ProcessCommandLine( const StringList& ) const override;
};

// ----------------------------------------------------------------------------

PCL_BEGIN_LOCAL
extern GHSProcess* TheGHSProcess;
PCL_END_LOCAL

// ----------------------------------------------------------------------------

} // pcl

#endif   // __GHSProcess_h

// ----------------------------------------------------------------------------
// EOF GHSProcess.h

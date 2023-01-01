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
// GHSModule.cpp
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

#define MODULE_VERSION_MAJOR     3
#define MODULE_VERSION_MINOR     0
#define MODULE_VERSION_REVISION  2
#define MODULE_VERSION_BUILD     0
#define MODULE_VERSION_LANGUAGE  eng
#define MODULE_VERSION_STATUS    release

#define MODULE_RELEASE_YEAR      2023
#define MODULE_RELEASE_MONTH     01
#define MODULE_RELEASE_DAY       01

#include "GHSModule.h"
#include "GHSProcess.h"
#include "GHSInterface.h"

namespace pcl
{

// ----------------------------------------------------------------------------

GHSModule::GHSModule()
{
}

// ----------------------------------------------------------------------------

const char* GHSModule::Version() const
{
   return PCL_MODULE_VERSION_S(MODULE_VERSION_MAJOR,
                               MODULE_VERSION_MINOR,
                               MODULE_VERSION_REVISION,
                               MODULE_VERSION_BUILD,
                               MODULE_VERSION_LANGUAGE,
                               MODULE_VERSION_STATUS);
}

// ----------------------------------------------------------------------------

IsoString GHSModule::Name() const
{
   /*
    * Replace with the actual name of this module. Must be unique.
    */
   return "GeneralizedHyperbolicStretch";
}

// ----------------------------------------------------------------------------

String GHSModule::Description() const
{
   return "PixInsight GHS Process Module";
}

// ----------------------------------------------------------------------------

String GHSModule::Company() const
{
   return "www.ghsastro.co.uk";
}

// ----------------------------------------------------------------------------

String GHSModule::Author() const
{
   return "Mike Cranfield";
}

// ----------------------------------------------------------------------------

String GHSModule::Copyright() const
{
   return "Copyright (c) 2022,2023 Mike Cranfield";
}

// ----------------------------------------------------------------------------

String GHSModule::TradeMarks() const
{
   return "PixInsight";
}

// ----------------------------------------------------------------------------

String GHSModule::OriginalFileName() const
{
#ifdef __PCL_FREEBSD
   return "GeneralizedHyperbolicStretch-pxm.so";
#endif
#ifdef __PCL_LINUX
   return "GeneralizedHyperbolicStretch-pxm.so";
#endif
#ifdef __PCL_MACOSX
   return "GeneralizedHyperbolicStretch-pxm.dylib";
#endif
#ifdef __PCL_WINDOWS
   return "GeneralizedHyperbolicStretch-pxm.dll";
#endif
}

// ----------------------------------------------------------------------------

void GHSModule::GetReleaseDate( int& year, int& month, int& day ) const
{
   year  = MODULE_RELEASE_YEAR;
   month = MODULE_RELEASE_MONTH;
   day   = MODULE_RELEASE_DAY;
}

// ----------------------------------------------------------------------------

} // pcl

// ----------------------------------------------------------------------------

/*
 * Module installation routine.
 *
 * If this routine is defined as a public symbol in a module, the PixInsight
 * core application will call it just after loading and initializing the module
 * shared object or dynamic-link library.
 *
 * The mode argument specifies the kind of installation being performed by the
 * core application. See the pcl::InstallMode namespace for more information.
 */
PCL_MODULE_EXPORT int InstallPixInsightModule( int mode )
{
   /*
    * When the PixInsight application installs this module, we just have to
    * instantiate the meta objects describing it.
    */
   new pcl::GHSModule;

   /*
    * The mode argument tells us which kind of installation is being requested
    * by the PixInsight application. Incomplete installation requests only need
    * module descriptions.
    */
   if ( mode == pcl::InstallMode::FullInstall )
   {
      new pcl::GHSProcess;
      new pcl::GHSInterface;
   }

   /*
    * Return zero to signal successful installation
    */
   return 0;
}

// ----------------------------------------------------------------------------
// EOF GHSModule.cpp

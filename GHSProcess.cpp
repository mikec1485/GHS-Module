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
// GHSProcess.cpp
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

#include "GHSInstance.h"
#include "GHSInterface.h"
#include "GHSParameters.h"
#include "GHSProcess.h"

#include <pcl/Arguments.h>
#include <pcl/Console.h>
#include <pcl/Exception.h>
#include <pcl/View.h>

namespace pcl
{

// ----------------------------------------------------------------------------

GHSProcess* TheGHSProcess = nullptr;

// ----------------------------------------------------------------------------

GHSProcess::GHSProcess()
{
   TheGHSProcess = this;

   /*
    * Instantiate process parameters.
    */
    new GHSST( this );
    new GHSSC( this );
    new GHSInv( this );
    new GHS( this );
    new GHSb( this );
    new GHSSP( this );
    new GHSHP( this );
    new GHSLP( this );
    new GHSBP( this );
    new GHSWP( this );
    new GHSCB( this );
    new GHSCT( this );
    new GHSRGBWS( this );
 
}

// ----------------------------------------------------------------------------

IsoString GHSProcess::Id() const
{
   return "GeneralizedHyperbolicStretch";
}

// ----------------------------------------------------------------------------

IsoString GHSProcess::Categories() const
{
   return IsoString("IntensityTransformations");
}

// ----------------------------------------------------------------------------

uint32 GHSProcess::Version() const
{
   return 0x100;
}

// ----------------------------------------------------------------------------

String GHSProcess::Description() const
{
   return
   "<html>"
   "<p>"
    "Applies intensity transformations to images using various mathematical forms adapted to allow  significant user control."
   "</p>"
   "</html>";
}

// ----------------------------------------------------------------------------

String GHSProcess::IconImageSVGFile() const
{
   return "@module_icons_dir/GHS.svg";
}

// ----------------------------------------------------------------------------

ProcessInterface* GHSProcess::DefaultInterface() const
{
   return TheGHSInterface;
}

// ----------------------------------------------------------------------------

ProcessImplementation* GHSProcess::Create() const
{
   return new GHSInstance( this );
}

// ----------------------------------------------------------------------------

ProcessImplementation* GHSProcess::Clone( const ProcessImplementation& p ) const
{
   const GHSInstance* instance = dynamic_cast<const GHSInstance*>( &p );
   return (instance != nullptr) ? new GHSInstance( *instance ) : nullptr;
}

// ----------------------------------------------------------------------------

bool GHSProcess::CanProcessCommandLines() const
{
   return true;
}

// ----------------------------------------------------------------------------

static void ShowHelp()
{
   Console().Write(
"<raw>"
"Usage: GeneralizedHyperbolicStretch [<arg_list>] [<view_list>]"
"\n"
"\n--interface"
"\n"
"\n      Launches the interface of this process."
"\n"
"\n--help"
"\n"
"\n      Displays this help and exits."
"</raw>" );
}

int GHSProcess::ProcessCommandLine( const StringList& argv ) const
{
   ArgumentList arguments = ExtractArguments( argv, ArgumentItemMode::AsViews, ArgumentOption::AllowWildcards );

   GHSInstance instance( this );

   bool launchInterface = false;
   int count = 0;

   for ( const Argument& arg : arguments )
   {
      if ( arg.IsNumeric() )
      {
         throw Error( "Unknown numeric argument: " + arg.Token() );
      }
      else if ( arg.IsString() )
      {
         throw Error( "Unknown string argument: " + arg.Token() );
      }
      else if ( arg.IsSwitch() )
      {
         throw Error( "Unknown switch argument: " + arg.Token() );
      }
      else if ( arg.IsLiteral() )
      {
         // These are standard parameters that all processes should provide.
         if ( arg.Id() == "-interface" )
            launchInterface = true;
         else if ( arg.Id() == "-help" )
         {
            ShowHelp();
            return 0;
         }
         else
            throw Error( "Unknown argument: " + arg.Token() );
      }
      else if ( arg.IsItemList() )
      {
         ++count;

         if ( arg.Items().IsEmpty() )
            throw Error( "No view(s) found: " + arg.Token() );

         for ( StringList::const_iterator j = arg.Items().Begin(); j != arg.Items().End(); ++j )
         {
            View v = View::ViewById( *j );
            if ( v.IsNull() )
               throw Error( "No such view: " + *j );
            instance.LaunchOn( v );
         }
      }
   }

   if ( launchInterface )
      instance.LaunchInterface();
   else if ( count == 0 )
   {
      if ( ImageWindow::ActiveWindow().IsNull() )
         throw Error( "There is no active image window." );
      instance.LaunchOnCurrentView();
   }

   return 0;
}

// ----------------------------------------------------------------------------

} // pcl

// ----------------------------------------------------------------------------
// EOF GHSProcess.cpp

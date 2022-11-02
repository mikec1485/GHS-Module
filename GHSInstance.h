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
// GHSInstance.h
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


#ifndef __GHSInstance_h
#define __GHSInstance_h

#include <pcl/MetaParameter.h> // pcl_bool, pcl_enum
#include <pcl/ProcessImplementation.h>
#include <pcl/Console.h>
#include <pcl/ReadoutOptions.h>

#include "GHSParameters.h"

namespace pcl
{

// ----------------------------------------------------------------------------

class GHSInstance : public ProcessImplementation
{
public:

    GHSInstance( const MetaProcess* );
    GHSInstance( const GHSInstance& );

    void Assign( const ProcessImplementation& ) override;
    UndoFlags UndoMode( const View& ) const override;
    bool CanExecuteOn( const View&, pcl::String& whyNot ) const override;
    bool CanExecuteOn( const ImageVariant&, pcl::String& whyNot ) const override;
    bool ExecuteOn( View& ) override;
    bool ExecuteOn( ImageVariant&, const IsoString& ) override;
    void* LockParameter( const MetaParameter*, size_type tableRow ) override;
    void Preview( UInt16Image& );
    void TransformHistogram( Histogram& , const Histogram& );

private:

    //-------------------------------------------------------
    
    void Transform( double& value ) const
    {
        if (!IsIdentityTransformation())   // no need to update value if identity stretch
        {            
         if (m_flags.Linear)   // Linear stretch
         {
             value = (value <= p_BP) ? 0.0 : ((value >= p_WP ) ? 1.0 : (value - p_BP)/(p_WP - p_BP));
             
             // A linear stretch is not truly invertible so we never allow this
             //if (m_flags.Inverse)
             //     value = value * (p_WP - p_BP) + p_BP;
             //else
             //     value = (value <= p_BP) ? 0.0 : ((value >= p_WP ) ? 1.0 : (value - p_BP)/(p_WP - p_BP));
         }
         else   // Non-linear stretch
         {
             if ((p_SP < p_LP) || (p_SP > p_HP))
                 return;
             
             if (m_flags.GHSLog)
             {
                 double LPT = p_LP;
                 double SPT = p_SP;
                 double HPT = p_HP;

                 if       (value < LPT)  {value = m_coeffs.a1 + m_coeffs.b1 * value;}
                 else if  (value < SPT)  {value = m_coeffs.a2 + m_coeffs.b2 * Ln(m_coeffs.c2 + m_coeffs.d2 * value);}
                 else if  (value < HPT)  {value = m_coeffs.a3 + m_coeffs.b3 * Ln(m_coeffs.c3 + m_coeffs.d3 * value);}
                 else                    {value = m_coeffs.a4 + m_coeffs.b4 * value;}
             }
             
             if (m_flags.GHSLogInv)
             {
                 double LPT = m_coeffs.a1 + m_coeffs.b1 * p_LP;
                 double SPT = m_coeffs.a2 + m_coeffs.b2 * Ln(m_coeffs.c2 + m_coeffs.d2 * p_SP);
                 double HPT = m_coeffs.a4 + m_coeffs.b4 * p_HP;

                 if       (value < LPT)  {value = (value - m_coeffs.a1) / m_coeffs.b1;}
                 else if  (value < SPT)  {value = (Exp((value - m_coeffs.a2) / m_coeffs.b2) - m_coeffs.c2) / m_coeffs.d2;}
                 else if  (value < HPT)  {value = (Exp((value - m_coeffs.a3) / m_coeffs.b3) - m_coeffs.c3) / m_coeffs.d3;}
                 else                    {value = (value - m_coeffs.a4) / m_coeffs.b4;}
             }
             
             if (m_flags.GHSExp)
             {
                 double LPT = p_LP;
                 double SPT = p_SP;
                 double HPT = p_HP;

                 if       (value < LPT)  {value = m_coeffs.a1 + m_coeffs.b1 * value;}
                 else if  (value < SPT)  {value = m_coeffs.a2 + m_coeffs.b2 * exp(m_coeffs.c2 + m_coeffs.d2 * value);}
                 else if  (value < HPT)  {value = m_coeffs.a3 + m_coeffs.b3 * exp(m_coeffs.c3 + m_coeffs.d3 * value);}
                 else                    {value = m_coeffs.a4 + m_coeffs.b4 * value;}
             }
             
             if (m_flags.GHSExpInv)
             {
                 double LPT = m_coeffs.a1 + m_coeffs.b1 * p_LP;
                 double SPT = m_coeffs.a2 + m_coeffs.b2 * Exp(m_coeffs.c2 + m_coeffs.d2 * p_SP);
                 double HPT = m_coeffs.a4 + m_coeffs.b4 * p_HP;

                 if       (value < LPT)  {value = (value - m_coeffs.a1) / m_coeffs.b1;}
                 else if  (value < SPT)  {value = (Ln((value - m_coeffs.a2) / m_coeffs.b2) - m_coeffs.c2) / m_coeffs.d2;}
                 else if  (value < HPT)  {value = (Ln((value - m_coeffs.a3) / m_coeffs.b3) - m_coeffs.c3) / m_coeffs.d3;}
                 else                    {value = (value - m_coeffs.a4) / m_coeffs.b4;}
             }
             
             if (m_flags.GHSHyp || m_flags.GHSInt)
             {
                 double LPT = p_LP;
                 double SPT = p_SP;
                 double HPT = p_HP;

                 if       (value < LPT)  {value = m_coeffs.a1 + m_coeffs.b1 * value;}
                 else if  (value < SPT)  {value = m_coeffs.a2 + m_coeffs.b2 * Pow((m_coeffs.c2 + m_coeffs.d2 * value), m_coeffs.e2);}
                 else if  (value < HPT)  {value = m_coeffs.a3 + m_coeffs.b3 * Pow((m_coeffs.c3 + m_coeffs.d3 * value), m_coeffs.e3);}
                 else                    {value = m_coeffs.a4 + m_coeffs.b4 * value;}
             }
             
             if (m_flags.GHSHypInv || m_flags.GHSIntInv)
             {
                 double LPT = m_coeffs.a1 + m_coeffs.b1 * p_LP;
                 double SPT = m_coeffs.a2 + m_coeffs.b2 * Pow((m_coeffs.c2 + m_coeffs.d2 * p_SP), m_coeffs.e2);
                 double HPT = m_coeffs.a4 + m_coeffs.b4 * p_HP;

                 if       (value < LPT)  {value = (value - m_coeffs.a1) / m_coeffs.b1;}
                 else if  (value < SPT)  {value = (Pow((value - m_coeffs.a2)/m_coeffs.b2, 1/m_coeffs.e2) - m_coeffs.c2) / m_coeffs.d2;}
                 else if  (value < HPT)  {value = (Pow((value - m_coeffs.a3)/m_coeffs.b3, 1/m_coeffs.e3) - m_coeffs.c3) / m_coeffs.d3;}
                 else                    {value = (value - m_coeffs.a4) / m_coeffs.b4;}
             }
                    
             if (m_flags.GHSMtf)
             {
                 double LPT = p_LP;
                 double SPT = p_SP;
                 double HPT = p_HP;

                 if       (value < LPT)  {value = m_coeffs.a1 + m_coeffs.b1 * value;}
                 else if  (value < SPT)  {value = m_coeffs.a2 + (m_coeffs.b2 * value + m_coeffs.c2)/(m_coeffs.d2 * value + m_coeffs.e2);}
                 else if  (value < HPT)  {value = m_coeffs.a3 + (m_coeffs.b3 * value + m_coeffs.c3)/(m_coeffs.d3 * value + m_coeffs.e3);}
                 else                    {value = m_coeffs.a4 + m_coeffs.b4 * value;}
             }
             
             if (m_flags.GHSMtfInv)
             {
                 double LPT = m_coeffs.a1 + m_coeffs.b1 * p_LP;
                 double SPT = m_coeffs.a2 + (m_coeffs.b2 * p_SP + m_coeffs.c2) / (m_coeffs.d2 * p_SP + m_coeffs.e2);
                 double HPT = m_coeffs.a4 + m_coeffs.b4 * p_HP;
                 
                 if       (value < LPT)   {value = (value - m_coeffs.a1) / m_coeffs.b1;}
                 else if  (value < SPT)   {value = (m_coeffs.c2 - m_coeffs.e2 * (value - m_coeffs.a2))/(m_coeffs.d2 * (value - m_coeffs.a2) - m_coeffs.b2);}
                 else if  (value < HPT)   {value = (m_coeffs.c3 - m_coeffs.e3 * (value - m_coeffs.a3))/(m_coeffs.d3 * (value - m_coeffs.a3) - m_coeffs.b3);}
                 else                 {value = (value - m_coeffs.a4) / m_coeffs.b4;}
             }
             
             if (m_flags.GHSAsh)
             {
                 double LPT = p_LP;
                 double SPT = p_SP;
                 double HPT = p_HP;

                 if       (value < LPT)  {value = m_coeffs.a1 + m_coeffs.b1 * value;}
                 else if  (value < SPT)  {value = m_coeffs.a2 + m_coeffs.b2 * Ln(m_coeffs.c2 * (value - m_coeffs.e2) + Sqrt(m_coeffs.d2 * (value - m_coeffs.e2) * (value - m_coeffs.e2) + 1));}
                 else if  (value < HPT)  {value = m_coeffs.a3 + m_coeffs.b3 * Ln(m_coeffs.c3 * (value - m_coeffs.e3) + Sqrt(m_coeffs.d3 * (value - m_coeffs.e3) * (value - m_coeffs.e3) + 1));}
                 else                    {value = m_coeffs.a4 + m_coeffs.b4 * value;}
             }
             
             if (m_flags.GHSAshInv)
             {
                 double LPT = m_coeffs.a1 + m_coeffs.b1 * p_LP;
                 double SPT = m_coeffs.a2 + m_coeffs.b2 * Ln(m_coeffs.c2 * (p_SP - m_coeffs.e2) + Sqrt(m_coeffs.d2 * (p_SP - m_coeffs.e2) * (p_SP - m_coeffs.e2) + 1));
                 double HPT = m_coeffs.a4 + m_coeffs.b4 * p_HP;

                 if       (value < LPT)  {value = (value - m_coeffs.a1) / m_coeffs.b1;}
                 else if  (value < SPT)  {double expVal = Exp((m_coeffs.a2 - value) / m_coeffs.b2); value = m_coeffs.e2 - (expVal - (1 / expVal)) / (2 * m_coeffs.c2);}
                 else if  (value < HPT)  {double expVal = Exp((m_coeffs.a3 - value) / m_coeffs.b3); value = m_coeffs.e3 - (expVal - (1 / expVal)) / (2 * m_coeffs.c3);}
                 else                    {value = (value - m_coeffs.a4) / m_coeffs.b4;}
             }
         }
        }
        value = Range(value, 0.0, 1.0);
    }
    
    //-------------------------------------------------------
    
    struct Flags
    {
        bool GHSType = false;
        bool GHSTypeInv = false;
        bool Linear = false;
        bool Inverse = false;
        bool GHSLog = false;
        bool GHSExp = false;
        bool GHSHyp = false;
        bool GHSInt = false;
        bool GHSMtf = false;
        bool GHSAsh = false;
        bool GHSLogInv = false;
        bool GHSExpInv = false;
        bool GHSHypInv = false;
        bool GHSIntInv = false;
        bool GHSMtfInv = false;
        bool GHSAshInv = false;

       /*!
        * Default constructor.
        */
       Flags() = default;

       /*!
        * Copy constructor.
        */
       Flags( const Flags& ) = default;

       /*!
        * Constructs a new %Flags object initialized for the parameters of the
        * specified GHSTransformation \a H.
        */
       Flags( GHSInstance& G )
       {
           Linear = (G.p_ST == GHSST::ST_Linear);
           Inverse = G.IsInvertible() && G.p_Inv;
           
           GHSType = (G.p_ST == GHSST::ST_GeneralisedHyperbolic) && (!Inverse);
           GHSLog = (GHSType) && (G.p_b == -1);
           GHSExp = (GHSType) && (G.p_b == 0);
           GHSInt = (GHSType) && (G.p_b < 0) && (!GHSLog);
           GHSHyp = (GHSType) && (G.p_b > 0);
           GHSMtf = (G.p_ST == GHSST::ST_MidtonesTransfer) && (!Inverse);
           GHSAsh = (G.p_ST == GHSST::ST_Arcsinh) && (!Inverse);
           
           GHSTypeInv = (G.p_ST == GHSST::ST_GeneralisedHyperbolic) && (Inverse);
           GHSLogInv = (GHSTypeInv) && (G.p_b == -1);
           GHSExpInv = (GHSTypeInv) && (G.p_b == 0);
           GHSIntInv = (GHSTypeInv) && (G.p_b < 0) && (!GHSLogInv);
           GHSHypInv = (GHSTypeInv) && (G.p_b > 0);
           GHSMtfInv = (G.p_ST == GHSST::ST_MidtonesTransfer) && (Inverse);
           GHSAshInv = (G.p_ST == GHSST::ST_Arcsinh) && (Inverse);
       }
    };

    /*!
     * Returns the set of flags characterizing this GHS transformation.
     */
    Flags TransformationFlags() const
    {
       return m_flags;
    }
    
    //-------------------------------------------------------

    struct Coeffs
    {
        double a1 = 0.0;
        double b1 = 0.0;

        double a2 = 0.0;
        double b2 = 0.0;
        double c2 = 0.0;
        double d2 = 0.0;
        double e2 = 0.0;

        double a3 = 0.0;
        double b3 = 0.0;
        double c3 = 0.0;
        double d3 = 0.0;
        double e3 = 0.0;

        double a4 = 0.0;
        double b4 = 0.0;
        
        /*!
         * Default constructor.
         */
        Coeffs() = default;

        /*!
         * Copy constructor.
         */
        Coeffs( const Coeffs& ) = default;

        /*!
         * Constructs a new %Coeffs object initialized for the parameters of the
         * specified GHSTransformation \a G.
         */
        Coeffs( GHSInstance& G )
        {
            double orgD = G.p_D;
            double D = Exp(orgD) - 1.0;
            double B = G.p_b;
            double SP = G.p_SP;
            double LP = Min(SP, G.p_LP);
            double HP = Max(SP, G.p_HP);
            
            //---------------------------
            // GHS Logarithmic - (b = -1)|
            //---------------------------
            if ( G.m_flags.GHSLog ||  G.m_flags.GHSLogInv )
            {
                double qlp = -1.0 * Ln(1.0 + D * (SP - LP));
                double q0 = qlp - D * LP / (1.0 + D * (SP - LP));
                double qwp = Ln(1.0 + D * (HP - SP));
                double q1 = qwp + D * (1.0 - HP) / (1.0 + D * (HP - SP));
                double q = 1.0 / (q1 - q0);

                // derive coefficients for x < LP
                a1 = 0.0;
                b1 = D / (1.0 + D * (SP - LP)) * q;

                // derive coefficients for x < SP
                a2 = (-q0) * q;
                b2 = -q ;
                c2 = 1.0 + D * SP;
                d2 = -D;
                e2 = 0.0;

                // derive coefficients for SP <= x <= HP
                a3 = (-q0) * q;
                b3 = q;
                c3 = 1.0 - D * SP;
                d3 = D;
                e3 = 0.0;

                // derive coefficients for x > HP
                a4 = (qwp - q0 - D * HP / (1.0 + D * (HP - SP))) * q;
                b4 = q * D / (1.0 + D * (HP - SP));
            }

            //-----------------------
            // GHS Integral - (b < 0)|
            //-----------------------
            if ( G.m_flags.GHSInt ||  G.m_flags.GHSIntInv )
            {
                double qlp = -(1.0 - Pow((1.0 - D * B * (SP - LP)), (B + 1.0) / B)) / (B + 1);
                double q0 = qlp - D * LP * (Pow((1.0 - D * B * (SP - LP)), 1.0 / B));
                double qwp = -(Pow((1.0 - D * B * (HP - SP)), (B + 1.0) / B) - 1.0) / (B + 1);
                double q1 = qwp + D * (1.0 - HP) * (Pow((1.0 - D * B * (HP - SP)), 1.0 / B));
                double q = 1.0 / (q1 - q0);

                // derive coefficients for x < LP
                a1 = 0.0;
                b1 = D * Pow(1.0 - D * B * (SP - LP), 1.0 / B) * q;

                // derive coefficients for LP <= x < SP
                a2 = -(1 / (B + 1) + q0) * q;
                b2 = q / ( B + 1);
                c2 = 1.0 - D * B * SP;
                d2 = D * B;
                e2 = (B + 1.0) / B;

                // derive coefficients for SP <= x <= HP
                a3 = (1 / (B + 1) - q0) * q;
                b3 = -q / (B + 1);
                c3 = 1.0 + D * B * SP;
                d3 = -D * B;
                e3 = (B + 1.0) / B;

                // derive coefficients for x > HP
                a4 = (qwp - q0 - D * HP * Pow((1.0 - D * B * (HP - SP)), 1.0 / B)) * q;
                b4 = D * Pow((1.0 - D * B * (HP - SP)), 1.0 / B) * q;
            }

            //--------------------------
            // GHS Exponential - (b = 0)|
            //--------------------------
            if ( G.m_flags.GHSExp ||  G.m_flags.GHSExpInv )
            {
                double qlp = Exp(-D * (SP - LP));
                double q0 = qlp - D * LP * qlp;
                double qwp = 2.0 - Exp(-D * (HP - SP));
                double q1 = qwp + D * (1.0 - HP) * (2.0 - qwp);
                double q = 1.0 / (q1 - q0);

                // derive coefficients for x < LP
                a1 = 0.0;
                b1 = D * qlp * q;

                // derive coefficients for LP <= x < SP
                a2 = -q0 * q;
                b2 = q;
                c2 = -D * SP;
                d2 = D;
                e2 = 0.0;

                // derive coefficients for SP <= x <= HP
                a3 = (2.0 - q0) * q;
                b3 = -q;
                c3 = D * SP;
                d3 = -D;
                e3 = 0.0;

                // derive coefficients for x > HP
                a4 = (qwp - q0 - D * HP * (2.0 - qwp)) * q;
                b4 = D * (2.0 - qwp) * q;
            }

            //----------------------------------
            // GHS Hyperbolic/Harmonic - (b > 0)|
            //----------------------------------
            if ( G.m_flags.GHSHyp ||  G.m_flags.GHSHypInv )
            {
                double qlp = Pow((1 + D * B * (SP - LP)), -1.0 / B);
                double q0 = qlp - D * LP * Pow((1 + D * B * (SP - LP)), -(1.0 + B) / B);
                double qwp = 2.0 - Pow(1.0 + D * B * (HP - SP), -1.0 / B);
                double q1 = qwp + D * (1.0 - HP) * Pow((1.0 + D * B * (HP - SP)), -(1.0 + B) / B);
                double q = 1.0 / (q1 - q0);

                // derive coefficients for x < LP
                a1 = 0.0;
                b1 = D * Pow((1 + D * B * (SP - LP)), -(1.0 + B) / B) * q;

                // derive coefficients for LP <= x < SP
                a2 = -q0 * q;
                b2 = q;
                c2 = 1.0 + D * B * SP;
                d2 = -D * B;
                e2 = -1.0 / B;

                // derive coefficients for SP <= x <= HP
                a3 = (2.0 - q0) * q;
                b3 = -q;
                c3 = 1.0 - D * B * SP;
                d3 = D * B;
                e3 = -1.0 / B;

                // derive coefficients for x > HP
                a4 = (qwp - q0 - D * HP * Pow((1.0 + D * B * (HP - SP)), -(B + 1.0) / B)) * q;
                b4 = (D * Pow((1.0 + D * B * (HP - SP)), -(B + 1.0) / B)) * q;
            }
            
            // MTF stretch (Histogram Transformation)
            if ( G.m_flags.GHSMtf || G.m_flags.GHSMtfInv )
            {
                double m = 1 / (2 * (D + 1));
                double zLPSP = (1 - 2 * m)*(LP - SP) - m;
                double qlp = (m - 1) * (LP - SP) / zLPSP;
                double q0 =  qlp + LP * (m - 1) * m / (zLPSP * zLPSP);
                double qwp = (m - 1) * (HP - SP) / ((2 * m - 1) * (HP - SP) - m);
                double zSPHP = (2 * m - 1)*(HP - SP) - m;
                double q1 =  qwp + (HP - 1) * (m - 1) * m / (zSPHP * zSPHP);
                double q = 1.0 / (q1 - q0);

                // derive coefficients for x < LP
                a1 = 0.0;
                b1 = m * (1 - m) * q / (zLPSP * zLPSP);

                // derive coefficients for LP <= x <= SP
                a2 = -q0 * q;
                b2 = (m-1) * q;
                c2 = b2 * (-SP);
                d2 = (1 - 2 * m);
                e2 = -d2 * SP - m;

                // derive coefficients for LP <= x <= SP
                a3 = -q0 * q;
                b3 = (m - 1) * q;
                c3 = b3 * (-SP);
                d3 = (2 * m - 1);
                e3 = -d3 * SP - m;

                // derive coefficients for x > HP
                a4 = (qwp - HP * (1 - m) * m / (zSPHP * zSPHP) - q0) * q;
                b4 = -m * (m - 1) * q / (zSPHP * zSPHP);
            }
            
            // Arcsinh Stretch
            if ( G.m_flags.GHSAsh || G.m_flags.GHSAshInv )
            {
                double powSPLP = Sqrt(D * D * (SP - LP) * (SP - LP) + 1);
                double powHPSP = Sqrt(D * D * (HP - SP) * (HP - SP) + 1);
                double qlp = - Ln(D * (SP - LP) + powSPLP);
                double q0 = qlp - LP * D / powSPLP;
                double qwp = Ln(D * (HP - SP) + powHPSP);
                double q1 = qwp + (1.0 - HP) * D / powHPSP;
                double q = 1.0 / (q1 - q0);

                // derive coefficients for x < LP
                a1 = 0.0;
                b1 = D * q / powSPLP;

                // derive coefficients for LP <= x < SP
                a2 = -q0 * q;
                b2 = -q;
                c2 = -D;
                d2 = D * D;
                e2 = SP;

                // derive coefficients for SP <= x <= HP
                a3 = -q0 * q;
                b3 = q;
                c3 = D;
                d3 = D * D;
                e3 = SP;

                // derive coefficients for x > HP
                a4 = (qwp - HP * D / powHPSP - q0) * q;
                b4 = D * q / powHPSP;
            }
            
            // Linear stretch (coefficients not needed)
            if ( G.m_flags.Linear )
            {
                // derive coefficients for x < LP
                a1 = 0.0;
                b1 = 0.0;

                // derive coefficients for LP <= x < SP
                a2 = 0.0;
                b2 = 0.0;
                c2 = 0.0;
                d2 = 0.0;
                e2 = SP;

                // derive coefficients for SP <= x <= HP
                a3 = 0.0;
                b3 = 0.0;
                c3 = 0.0;
                d3 = 0.0;
                e3 = 0.0;

                // derive coefficients for x > HP
                a4 = 0.0;
                b4 = 0.0;
            }
        }
    };
    
    /*!
     * Returns the set of coefficients characterizing this GHS transformation.
     */
    Coeffs TransformationCoeffs() const
    {
       return m_coeffs;
    }
    
    //-------------------------------------------------------
    
    void UpdateFlags()
    {
        m_flags = Flags( *this );
    }
     
    void UpdateCoeffs()
    {
        m_coeffs = Coeffs( *this );
    }
    
    //-------------------------------------------------------
    
    pcl_enum    p_ST;
    pcl_enum    p_SC;
    pcl_bool   p_Inv;
    double      p_D;
    double      p_b;
    double      p_SP;
    double      p_LP;
    double      p_HP;
    double      p_BP;
    double      p_WP;
    double      p_CB;
    pcl_enum    p_CT;
    pcl_bool    p_RGBWS;
    
    
    
    Flags       m_flags;
    Coeffs      m_coeffs;
    
    //-------------------------------------------------------
    
    /*! #
     */
    void SetStretchType( int ST, bool recalc = true )
    {
        p_ST = pcl::Range( ST, 0, GHSST::ST_NumberOfItems - 1 );
        
        if (recalc)
        {
            UpdateFlags();
            UpdateCoeffs();
        }
       
    }

     /*! #
      */
     void SetStretchChannel( int SC, bool recalc = true )
     {
         p_SC = pcl::Range( SC, 0, GHSSC::SC_NumberOfItems - 1 );
         
         if (recalc)
         {
             //UpdateFlags();
             //UpdateCoeffs();
         }
     }

     void SetInverse( bool Inv, bool recalc = true )
     {
         p_Inv = Inv;
         
         if (recalc)
         {
             UpdateFlags();
             //UpdateCoeffs();
         }
     }

     /*! #
     */
    void SetStretchFactor( double D, bool recalc = true )
    {
        p_D = pcl::Max(0.0, D);
        
        if (recalc)
        {
            //UpdateFlags();
            UpdateCoeffs();
        }
    }

    /*! #
    */
   void SetLocalIntensity( double b, bool recalc = true )
   {
       p_b = b;
       
       if (recalc)
       {
           UpdateFlags();
           UpdateCoeffs();
       }
   }

   /*! #
     */
    void SetSymmetryPoint( double SP, bool recalc = true )
    {
        p_LP = pcl::Max(0.0, pcl::Min(p_LP, SP));
        p_HP = pcl::Min(1.0, pcl::Max(SP, p_HP));
        p_SP = pcl::Range( SP, 0.0, 1.0 );
        
        if (recalc)
        {
            //UpdateFlags();
            UpdateCoeffs();
        }
    }

    /*! #
     */
    void SetShadowProtect( double LP, bool recalc = true )
    {
        p_LP = pcl::Range( LP, 0.0, p_SP );
        
        if (recalc)
        {
            //UpdateFlags();
            UpdateCoeffs();
        }
    }

     /*! #
      */
     void SetHighlightProtect( double HP, bool recalc = true )
     {
         p_HP = pcl::Range( HP, p_SP, 1.0 );
         
         if (recalc)
         {
             //UpdateFlags();
             UpdateCoeffs();
         }
     }

     /*! #
      */
     void SetBlackpoint( double BP, bool recalc = true )
     {
         double x = double(TheGHSWPParameter->Precision());
         p_BP = pcl::Min( BP, 1.0 );
         p_WP = pcl::Max( p_WP, p_BP + Pow10(-x));
         
         if (recalc)
         {
             //UpdateFlags();
             //UpdateCoeffs();
         }
     }

     /*! #
      */
     void SetWhitepoint( double WP, bool recalc = true )
     {
         double x = double(TheGHSBPParameter->Precision());
         p_WP = pcl::Max( WP, 0.0 );
         p_BP = pcl::Min( p_WP - Pow10(-x), p_BP );
         
         if (recalc)
         {
             //UpdateFlags();
             //UpdateCoeffs();
         }
     }

     /*! #
     */
    void SetColourBlend( double CB, bool recalc = true )
    {
        p_CB = pcl::Range( CB, 0.0, 1.0 );
        
        if (recalc)
        {
            //UpdateFlags();
            //UpdateCoeffs();
        }
    }

     /*! #
     */
    void SetClipType( int CT, bool recalc = true )
    {
        p_CT = pcl::Range( CT, 0, GHSCT::CT_NumberOfItems - 1 );
        
        if (recalc)
        {
            //UpdateFlags();
            //UpdateCoeffs();
        }
    }

    void SetRGBWS( bool rgbws, bool recalc = true )
    {
        p_RGBWS = rgbws;
        
        if (recalc)
        {
            //UpdateFlags();
            //UpdateCoeffs();
        }
    }

    bool IsIdentityTransformation() const
    {
       return ((p_ST < 3) && (p_D == 0.0)) ||
        ((p_ST == 3) && (p_BP == 0.0) && (p_WP == 1.0));
    }
    
    bool IsInvertible() const
    {
        bool saturation = (p_ST != GHSST::ST_Linear) && (p_SC == GHSSC::SC_Saturation);
        bool colour = (p_ST != GHSST::ST_Linear) && (p_SC == GHSSC::SC_Colour);
        bool linear = p_ST == GHSST::ST_Linear;
        return !(saturation || colour || linear);
    }
    
// ----------------------------------------------------------------------------


    friend class GHSEngine;
    friend class GHSProcess;
    friend class GHSInterface;
    friend struct Flags;
    friend struct Coeffs;
};

// ----------------------------------------------------------------------------


} // pcl

#endif   // __GHSInstance_h

// ----------------------------------------------------------------------------
// EOF GHSInstance.h

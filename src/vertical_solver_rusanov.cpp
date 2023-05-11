

#include "aether.h"
#include <cassert>
#include <math.h>

//This function modifies new_log_rho, new_log_NS, new_log_INS, new_vel_GD, new_temp, new_vert_vel
void advance_vertical_1stage(const arma_vec &log_rho, const arma_mat &log_NS, const arma_mat
                             &log_INS, const arma_mat &vel_GD, const arma_mat &i_vel, const arma_vec
                             &temp, const arma_mat &vert_vel, arma_vec &new_log_rho, arma_mat
                             &new_log_NS, arma_mat &new_log_INS, arma_mat &new_vel_GD, arma_vec
                             &new_temp, arma_mat &new_vert_vel) {
    
    //TODO: THESE VALUES ARE FROM OTHER FILES. DELETE LATER
    int iEast, iNorth, iUp;
    arma_vec InvRadialDistance_C;
    arma_vec mass; //TODO: Figure out where the mass vector comes from
    arma_vec dalt_f;
    precision_t Dt; //from MODGitm
    int nIonsAdvect; //in modPlanet (where is modPlanet)
    bool UseDamping;
    arma_vec vert_tau;
    precision_t Boltzmanns_Constant;
    bool UseCoriolis;
    precision_t MaximumVerticalVelocity;

    
    //TODO: Does fill::none need arma:: in front of it? I don't need it in constructor - delete later.
    arma_mat NS(nAlts+2, nSpecies, fill::none);
    arma_mat NS_small(nAlts, nSpecies, fill::none);
    arma_vec rho(nAlts+2, fill::none);
    arma_vec log_num(nAlts+2, fill::none);
    
    arma_vec grad_log_rho(nAlts, fill::none);
    arma_vec div_vel(nAlts, fill::none);
    arma_vec grad_temp(nAlts, fill::none);
    arma_vec grad_temp_kom(nAlts, fill::none);
    arma_vec diff_log_rho(nAlts, fill::none);
    arma_vec diff_temp(nAlts, fill::none);
    arma_vec grad_tmp(nAlts, fill::none); //NOTE: Different than grad_temp
    arma_vec diff_tmp(nAlts, fill::none);
    arma_vec diff_log_num(nAlts, fill::none);
    arma_vec grad_log_num(nAlts, fill::none);
    arma_vec divi_vel(nAlts, fill::none); //NOTE: Different than divi_vel
    
    arma_mat grad_vel_cd(nAlts, 3, fill::none);
    arma_mat diff_vel_cd(nAlts, 3, fill::none);
    arma_mat grad_ivel_cd(nAlts, 3, fill::none);
    arma_mat diff_ivel_cd(nAlts, 3, fill::none);
    
    arma_mat grad_log_NS(nAlts, nSpecies, fill::none);
    arma_mat diff_log_ns(nAlts, nSpecies, fill::none);
    arma_mat grad_vert_vel(nAlts, nSpecies, fill::none);
    arma_mat diff_vert_vel(nAlts, nSpecies, fill::none);
    arma_mat div_vert_vel(nAlts, nSpecies, fill::none);
    
    arma_mat grad_log_INS(nAlts, nIons, fill::none);
    arma_mat diff_log_INS(nAlts, nIons, fill::none);
    precision_t new_sum_rho;
    precision_t new_log_sum_rho;
    
    int i_alt, i_species, j_species, i_dim;
    
    arma_vec NT(nAlts+2, fill::none);
    
    arma_mat nvel(nAlts, nSpecies, fill::none);
    int n_filter, i_filter;
    double low_filter;
    
    //Parameters used for the sponge:
    //This sponge is used to dampen out spurious modes oscillating
    //between the bottom and top of the model
    
    int n_alts_sponge = 12;
    precision_t ksp, nu_sp, amp_sp;
    
    //Eddy diffusion variables:
    arma_mat grad_log_con_s(nAlts, nSpecies, fill:none);
    arma_mat log_con_s(nAlts+2, nSpecies, fill:none);
    arma_mat con_s(nAlts+2, nSpecies, fill:none);
    
    
    
    NS = arma::exp(log_NS);
    rho = arma::exp(log_rho);
    
    //find the sum of every row in NS, take its natural log, and add it to log_num:
    for (uword i = 0; i < log_num.n_rows; i++) {
        precision_t row_sum = 0;
        for (uword j = 0; j < NS.n_cols; j++) {
            row_sum += NS(i, j);
        }
        log_num(i) = log(row_sum);
    }
    
    n_filter = 10;
    NT = arma::exp(log_num);
    
    calc_rusanov_alts_rusanov(log_rho, grad_log_rho, diff_log_rho);
    calc_rusanov_alts_rusanov(log_num, grad_log_num, diff_log_num);
    calc_rusanov_alts_rusanov(temp, grad_temp, diff_temp);
    
    for (size_t i = 0; i < 3; i++) {
        calc_rusanov_alts_rusanov(vel_gd.col(i), grad_vel_cd.col(i), diff_vel_cd.col(i));
        calc_rusanov_alts_rusanov(i_vel.col(i), grad_ivel_cd.col(i), diff_ivel_cd.col(i));
    }
    
    //Add geometrical correction to gradient and obtain divergence:
    div_vel = grad_vel_cd.col(iUp) + 2 * vel_gd.col(iUp) * InvRadialDistance_C;
    div_vel = grad_ivel_cd.col(iUp) + 2 * i_vel.col(iUp) * InvRadialDistance_C;
    
    for (uword i = 0; i < nSpecies; i++) {
        calc_rusanov_alts_rusanov(log_NS.col(i), grad_tmp, diff_tmp);
        assert(arma::size(grad_log_NS.col(i)) == arma::size(grad_tmp));
        grad_log_NS.col(i) = grad_tmp;
        assert(arma::size(diff_log_NS.col(i)) == arma::size(diff_tmp));
        diff_log_NS.col(i) = diff_tmp;
        
        calc_rusanov_alts_rusanov(vert_vel.col(i), grad_tmp, diff_tmp);
        assert(arma::size(grad_vert_vel.col(i)) == arma::size(grad_tmp));
        grad_vert_vel.col(i) = grad_tmp;
        assert(arma::size(diff_vert_vel.col(i)) == arma::size(diff_tmp));
        diff_vert_vel.col(i) = diff_tmp;
        
        div_vert_vel.col(i) = grad_vert_vel.col(i) + 2 * vert_vel.col(i) * InvRadialDistance_C;
        
    }
    
    for (uword i = 0; i < nIons - 1; i++) {
        calc_rusanov_alts_rusanov(log_INS.col(i), grad_tmp, diff_tmp);
        grad_log_INS.col(i) = grad_tmp;
        diff_log_INS.col(i) = diff_tmp
    }
    //4th Order Gradients on a Non-Uniform Mesh (5-point Stencil)
    //Used for calculating the d(ln[Chi])/dr -> Log of the concentration gradient
    precision_t h1, h2, h3, h4;
    precision_t mesh_h2, mesh_h3, mesh_h4;
    precision_t mesh_coef_0, mesh_coef_1, mesh_coef_2, mesh_coef_3, mesh_coef_4;

    //Using Colegrove Method
    //Step 1: Calculator Ln(rho_{s}/Rho)
    for (uword i = 0; i < nSpecies; i++) {
        log_con_s.col(i) = log(mass(i) * NS.col(i) / rho);
    }
                               
    for (uword i = 0; i < nAlts - ; i++) {
        //4th order concentration gradient
        //on non-uniform mesh requires 5-pt stencil
        //TODO: Still need to figure this out
        h1 = dAlt_F(i - 1);
        h2 = dAlt_F(i);
        h3 = dAlt_F(i + 1);
        h4 = dAlt_F(i + 2);

        mesh_h2 = h2 + h1;
        mesh_h3 = h3 + h2 + h1;
        mesh_h4 = h4 + h3 + h2 + h1;

        mesh_coef_0 = (h2*h3*(h3+h4))/(h1*mesh_h2*mesh_h3*mesh_h4);
        mesh_coef_1 = -1.0*(mesh_h2*h3*(h3 + h4))/(h1*h2*(h2+h3)*(h2+h3+h4));
        mesh_coef_2 = (h2*h3*(h3+h4) + mesh_h2*h3*(h3+h4) -mesh_h2*h2*(h3+h4) -
                  mesh_h2*h2*h3)/(mesh_h2*h2*h3*(h3+h4));
        mesh_coef_3 = mesh_h2*h2*(h4 + h3)/(mesh_h3*(h2+h3)*h3*h4);
        mesh_coef_4 = -1.0*mesh_h2*h2*h3/(mesh_h4*(h2+h3+h4)*(h3+h4)*h4);

        for (uword ispecies = 0; i < nSpecies; i++) {
            grad_log_con_s(i_alt, iSpecies) = mesh_coef_0 * log_con_s(i_alt-2, ispecies)
            + mesh_coef_1 * log_con_s(i_alt-1, ispecies)
            + mesh_coef_2 * log_con_s(i_alt, ispecies)
            + mesh_coef_3 * log_con_s(i_alt+1, ispecies)
            + mesh_coef_4 * log_con_s(i_alt+2, ispecies);
        } //for inner  
    } //for outer

    
    amp_sp = 1 / (10 * Dt);
    ksp = n_alts_sponge + 1;
                               
    for (uword i_alt = 0; i_alt < nAlts; i_alt++) {
        for (uword ispecies = 0; ispecies < nSpecies; ispecies++) {
            new_log_NS(i_alt, ispecies) = log_NS(i_alt, ispecies) - Dt * (div_vert_vel(i_alt, ispecies)
            + vert_vel(i_alt, ispecies) * grad_log_NS(i_alt, ispecies) ) + Dt * diff_log_NS(i_alt, ispecies);
        }
        for (uword ispecies = 0; ispecies < nIonsAdvect; ispecies++) {
            new_log_INS(i_alt, ispecies) = log_INS(i_alt, ispecies) - Dt * (ivel(i_alt, ispecies)
            + vert_vel(i_alt, ispecies) * grad_log_INS(i_alt, ispecies) ) + Dt * diff_log_INS(i_alt, ispecies);
        }

        new_vel_GD(i_alt, iUp_) = 0;
        
        nu_sp = i_alt >= (nAlts - n_alts_sponge) ? amp_sp * (1 - arma::cos(M_PI * (ksp - (nAlts - i_alt))/ksp)) : 0;

        if (UseDamping) {
            vert_tau(i_alt) = 15 - 5 * (1 - arma::exp(-1 * Altitude_G(i_alt)/40000));
        }

        for (uword ispecies = 0; ispecies < nspecies) {

            new_vert_vel(i_alt, ispecies) = vert_vel(i_alt, ispecies) - 
                Dt * vert_vel(i_alt, iSpecies) * grad_vert_vel(i_alt, iSpecies) - 
                (vel_GD(i_alt, iNorth) * vel_GD(i_alt, iNorth) +
                vel_GD(i_alt, iEast) * vel_GD(i_alt, iEast)) *
                InvRadialDistance_C(i_alt) +
                temp(i_alt) * grad_log_NS(i_alt, ispecies) * Boltzmanns_Constant / mass(ispecies) +
                grad_temp(i_alt) * Boltzmanns_Constant / mass(ispecies) -
                Gravity_G(i_alt) + Dt * diff_vert_vel(i_alt, ispecies) -
                vert_vel(i_alt, ispecies) / vert_tau(i_alt);

            if (UseCoriolis) {
                new_vert_vel(i_alt, ispecies) = new_vert_vel(i_alt) +
                    Dt * (Centrifugal / InvRadialDistance_C(iAlt) Coriolis * vel_GD(i_alt, iEast));
            }

            //Thermal Diffusion Effects (For Light Species H2, H, and He) 
            // ThermalDiffCoefS is set in calc_rates
            // Note:  ThermalDiffCoefS is < 0.0 for light species
            // This forces light species into hot zones and heavy species into cold zones
            new_vert_vel(i_alt, ispecies) = new_vert_vel(i_alt, ispecies) - 
                Dt * (ThermalDiffCoefS(ispecies) * Boltzmanns_Constant * grad_temp(i_alt) / mass(ispecies));
        } //inner
    } //outer

    for (uword ialt = 0; ialt < nAlts; ialt++) {
        for (uword ispecies = 0; ispecies < nAlts; ispecies++) {
            new_vert_vel(ialt, ispecies) = std::min(MaximumVerticalVelocity, new_vert_vel(ialt, ispecies));
                new_vel_GD(ialt, iUp) = new_vel_GD(ialt, iUp) + new_vert_vel(ialt, ispecies) *
                (Mass(ispecies) * NS(ialt, ispecies) / rho(ialt));

        }
    }

    for (uword ialt = 0; ialt < nAlts; ialt++) {
        //dVphi/dt = - (V grad V)_phi
        new_vel_GD(ialt, iEast) = new_vel_GD(ialt, iEast) - 
            Dt * vel_GD(ialt, iUp) * grad_vel_cd(ialt, iEast) +
            Dt * diff_vel_cd(ialt, iEast);

        //dVtheta/dt = - (V grad V)_theta
        new_vel_GD(ialt, iNorth) = new_vel_GD(ialt, iNorth) - 
            Dt * vel_GD(ialt, iUp) * grad_vel_cd(ialt, iNorth) +
            Dt * diff_vel_cd(ialt, iNorth);

        new_temp(ialt) = new_temp(ialt) - 
        Dt * (vel_GD(ialt, iUp) * grad_temp(ialt) + (Gamma_1d(ialt) - 1) *
        (temp(ialt) * div_vel(ialt))) + Dt * diff_temp(ialt);
        
    }

    for (uword ialt = 0; ialt < nAlts; ialt++) {
        new_sum_rho = sum(Mass * exp(new_log_NS(ialt)));
        new_log_rho(ialt) = log(new_sum_rho);
    }




    
    

/*
 use ModGITM, only: &
        Dt, iEast_, iNorth_, iUp_, ThermalDiffCoefS
   use ModPlanet
   use ModSizeGitm
   use ModVertical, only : &
        EddyCoef_1D, Centrifugal, Coriolis, &
        MeanMajorMass_1d, Gamma_1d, InvRadialDistance_C, &
        KappaTemp_1d, &
        Gravity_G, Altitude_G,Cv_1D, dAlt_F
   use ModTime
   use ModInputs
   use ModConstants
 */
    
    
    
} //advance_vertical_1stage














                               



//modifies grad_var and diff_var
void calc_rusanov_alts_rusanov(const arma_vec &var, arma_vec &grad_var,const arma_vec &diff_var) {

    assert(var.size() == nAlts + 4); //TODO get nAlts in this function. nAlts member variable maybe?
    assert(diffVar.size() == nAlts);
    assert(gradVar.size() == nAlts);
    
    arma_vec varLeft(nAlts+1, 0);
    arma_vec varRight(nAlts+1, 0);
    arma_vec diffFlux(nAlts+1, 0);
    
    calc_facevalues_alts_rusanov(var, varLeft, varRight); //filling out varLeft and varRight
    
    //gradient based on averaged Left/Right values
    for (size_t i = 0; i < nAlts; i++) {
        gradVar[i] = (0.5 * varLeft[i+1] + varRight[i+1] - varLeft[i] - varRight[i])/dAlt_C[i]; //TODO get dAlt_C from modVertical
    } //for
    
    assert(cMax.size() > diffFlux.size()); //see comment in loop for explanation
    
    //Rusanov/Lax-Friedrichs diffusive term
    for (size_t i = 0; i < nAlts + 1; i++) {
        diffFlux[i] = 0.5 * std::max(cMax[i], cMax[i+1]) * (varRight[i] - varLeft[i]);
        //TODO get cMax array into this function. Also if cMax's size is the same or less than diffFlux the line above will cause segfault
    }
    
    for (size_t i = 0; i < nAlts; i++) {
        diffVar[i] = (diffFlux[i+1] - diffFlux[i]) / dAlt_C[i];
    }
    
    
} //calc_rusanov_alts_rusanov

void calc_facevalues_alts_rusanov(const arma_vec &var, arma_vec &varLeft, arma_vec &varRight) {

    assert(var.size() == nAlts + 4);
    assert(diffVar.size() == nAlts+1);
    assert(gradVar.size() == nAlts+1);
    
    precision_t dVarUp, dVarDown;
    arma_vec dVarLimited(nAlts + 2, 0);
    
    const precision_t factor1 = 0.6250000;
    const precision_t factor2 = 0.0416667;
    
    precision_t h;
    int i;
    
    //invDAlt_F is a arma_vec of size nAlts+3
    for (size_t i = 0; i < nAlts; i++) {
    
        h = invDAlt_F[i+1]*2;
        dVarUp = h * (factor1*(var[i-1]-var[i]) - factor2*(var[i+2]-Var[i-1]))
        
    }
    
} //calc_facevalues_alts_resanov

// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "../include/aether.h"


// -----------------------------------------------------------------------------
// Call (fortran) ionospheric electrodynamics if we enable this 
// in the cmake file
// -----------------------------------------------------------------------------



// -----------------------------------------------------------------------------
// Initialize the electrodynamics by reading the electrodynamics
// netCDF file.
// -----------------------------------------------------------------------------

Electrodynamics::Electrodynamics(Times time) {

  IsOk = true;

  bool isCompiled, isSet;
  int iError = 0;

  HaveElectrodynamicsFile = false;
  HaveFortranIe = false;
  isSet = false;

  #ifdef FORTRAN
    // Initialize IE components (reading in the data files):
    std::string ieDir = input.get_electrodynamics_dir();
    int* ieDir_iArray = copy_string_to_int(ieDir);
    std::string efield = input.get_potential_model();
    std::string aurora = input.get_diffuse_auroral_model();

    if (efield.length() == 0 & aurora.length() == 0) {
      HaveFortranIe = false;
    } else {
      int* efield_iArray = copy_string_to_int(efield);
      int* aurora_iArray = copy_string_to_int(aurora);
      ie_init_library(ieDir_iArray, efield_iArray, aurora_iArray, &iError);
      if (iError > 0) {
        report.print(0,"Error in setting fortran IE!");
        IsOk = false;
      } else {
        HaveFortranIe = true;
        isSet = true;
      }
    }
    isCompiled = true;
  #else
    isCompiled = false;
  #endif

  // If we don't set IE through Fortran, then try to set it through a file:

  if (!HaveFortranIe) {
    std::string electrodynamics_file = input.get_electrodynamics_file();
    if (electrodynamics_file.length() > 0 &
        electrodynamics_file != "none") {
      // This function sets HaveElectrodynamicsFile = true.
      read_netcdf_electrodynamics_file(electrodynamics_file);

      bool times_are_aligned = check_times(time.get_current(), time.get_end());

      if (!times_are_aligned) {
        IsOk = false;

        if (iProc == 0) {
          std::cout << "Times don't align with electrodynamics file! ";
          std::cout << "Please check this!\n";
        }
      }
    }
  }

  IsOk = sync_across_all_procs(IsOk);
}

// -----------------------------------------------------------------------------
// Update a bunch of indices and sent them to Fortran code.
// -----------------------------------------------------------------------------

void Electrodynamics::set_all_indices_for_ie(Times time,
                                              Indices &indices) {
  std::string function = "Electrodynamics::set_all_indices_for_ie";
  static int iFunction = -1;
  report.enter(function, iFunction);

  double time_now = time.get_current();
  int64_t iBy_ = indices.lookup_index_id("imfby");
  int64_t iBz_ = indices.lookup_index_id("imfbz");
  int64_t iVx_ = indices.lookup_index_id("swvx");
  int64_t iN_ = indices.lookup_index_id("swn");
  int64_t iAE_ = indices.lookup_index_id("ae");
  int64_t iAU_ = indices.lookup_index_id("au");
  int64_t iAL_ = indices.lookup_index_id("al");

  float imfby = indices.get_index(time_now, iBy_);
  float imfbz = indices.get_index(time_now, iBz_);
  float swv = indices.get_index(time_now, iVx_);
  float swn = indices.get_index(time_now, iN_);
  float ae = indices.get_index(time_now, iAE_);
  float au = indices.get_index(time_now, iAU_);
  float al = indices.get_index(time_now, iAL_);

  if (report.test_verbose(3)) {
    std::cout << "imf by : " << iBz_ << " " << imfby << "\n";
    std::cout << "imf bz : " << iBz_ << " " << imfbz << "\n";
    std::cout << "sw v : " << iVx_ << " " << swv << "\n";
    std::cout << "sw n : " << iN_ << " " << swn << "\n";
  }
  #ifdef FORTRAN
    ie_set_time(&time_now);
    ie_set_imfby(&imfby);
    ie_set_imfbz(&imfbz);
    ie_set_swv(&swv);
    ie_set_swn(&swn);
    ie_set_ae(&ae);
    ie_set_au(&au);
    ie_set_al(&al);
    ie_set_hp_from_ae(&ae);
  #endif

  report.exit(function);
  return;
}


// -----------------------------------------------------------------------------
// Update Electrodynamics
// -----------------------------------------------------------------------------

bool Electrodynamics::update(Planets planet,
                             Grid gGrid,
                             Times time,
                             Indices &indices,
                             Ions &ions) {


  std::string function = "Electrodynamics::update";
  static int iFunction = -1;
  report.enter(function, iFunction);

  if (HaveElectrodynamicsFile  || HaveFortranIe) {
    set_time(time.get_current());
    gGrid.calc_sza(planet, time);
    gGrid.calc_gse(planet, time);
    gGrid.calc_mlt();

    // Default is to set everything to zero:
    ions.potential_scgc.zeros();
    ions.eflux.zeros();
    ions.avee.ones();

    #ifdef FORTRAN
      if (HaveFortranIe) {
        report.print(3, "Using Fortran Electrodynamics!");
        set_all_indices_for_ie(time, indices);

        if (!IsAllocated) {
          int nXs = gGrid.get_nX();
          ie_set_nxs(&nXs);
          int nYs = gGrid.get_nY();
          ie_set_nys(&nYs);
          int64_t iTotal = nXs * nYs;
          mlt2d = static_cast<float*>(malloc(iTotal * sizeof(float)));
          lat2d = static_cast<float*>(malloc(iTotal * sizeof(float)));
          pot2d = static_cast<float*>(malloc(iTotal * sizeof(float)));
          eflux2d = static_cast<float*>(malloc(iTotal * sizeof(float)));
          avee2d = static_cast<float*>(malloc(iTotal * sizeof(float)));
          IsAllocated = true;
        }

        int64_t nZs = gGrid.get_nZ();
        int64_t iZ;

        int iError;

        for (iZ = 0; iZ < nZs; iZ++) {
          copy_mat_to_array(gGrid.magLocalTime_scgc.slice(iZ), mlt2d, true);
          copy_mat_to_array(gGrid.magLat_scgc.slice(iZ), lat2d, true);

          ie_set_mlts(mlt2d, &iError);
          ie_set_lats(lat2d, &iError);
          ie_update_grid(&iError);

          ie_get_potential(pot2d, &iError);
          copy_array_to_mat(pot2d, ions.potential_scgc.slice(iZ), true);

          if (iZ == nZs-1) {
            ie_get_electron_diffuse_aurora(eflux2d, avee2d, &iError);
            copy_array_to_mat(avee2d, ions.avee, true);
            copy_array_to_mat(eflux2d, ions.eflux, true);
          }
        }
      }
    #endif

    if (HaveElectrodynamicsFile) {
      report.print(3, "Setting electrodynamics from file!");
      auto electrodynamics_values =
        get_electrodynamics(gGrid.magLat_scgc,
                             gGrid.magLocalTime_scgc);
      ions.potential_scgc = std::get<0>(electrodynamics_values);
      ions.eflux = std::get<1>(electrodynamics_values);
      ions.avee = std::get<2>(electrodynamics_values);
    }
  } 

  report.exit(function);
  return true;
}

// -----------------------------------------------------------------------------
// Gets potential with generic interpolation scheme
// -----------------------------------------------------------------------------

arma_cube Electrodynamics::get_potential(arma_cube magLat,
                                         arma_cube magLocalTime) {
  arma_cube pot(magLat.n_rows, magLat.n_cols, magLat.n_slices);
  pot.zeros();

  if (HaveElectrodynamicsFile) {
    int time_pos = static_cast<int>(time_index);
    arma_mat e_potentials = input_electrodynamics[0].potential[time_pos];

    for (int i = 0; i < magLat.n_slices; ++i) {
      set_grid(magLat.slice(i) * cRtoD, magLocalTime.slice(i));
      pot.slice(i) = get_values(e_potentials, magLat.n_rows, magLat.n_cols);
    }
  }

  return pot;
}

// -----------------------------------------------------------------------------
// Gets energy flux with generic interpolation scheme
// -----------------------------------------------------------------------------

arma_mat Electrodynamics::get_eflux(arma_cube magLat,
                                    arma_cube magLocalTime) {

  arma_mat eflux(magLat.n_rows, magLat.n_cols);
  eflux.zeros();

  if (HaveElectrodynamicsFile) {
    int i = magLat.n_slices - 1;
    set_grid(magLat.slice(i) * cRtoD, magLocalTime.slice(i));
    int time_pos = static_cast<int>(time_index);
    arma_mat e_e_flux = input_electrodynamics[0].energy_flux[time_pos];
    eflux = get_values(e_e_flux, magLat.n_rows, magLat.n_cols);
  }

  return eflux;
}

// -----------------------------------------------------------------------------
// Gets average energy with generic interpolation scheme
// -----------------------------------------------------------------------------

arma_mat Electrodynamics::get_avee(arma_cube magLat,
                                   arma_cube magLocalTime) {
  arma_mat avee(magLat.n_rows, magLat.n_cols);
  avee.zeros();

  if (HaveElectrodynamicsFile) {
    int i = magLat.n_slices - 1;
    set_grid(magLat.slice(i) * cRtoD, magLocalTime.slice(i));
    int time_pos = static_cast<int>(time_index);
    arma_mat e_avee = input_electrodynamics[0].average_energy[time_pos];
    avee = get_values(e_avee, magLat.n_rows, magLat.n_cols);
  }

  return avee;
}

// -----------------------------------------------------------------------------
// Generic function to conduct 2d linear interpolation
// -----------------------------------------------------------------------------

arma_mat Electrodynamics::get_values(arma_mat matToInterpolateOn,
                                     int rows, int cols) {
  arma_mat slice(rows, cols);
  slice.zeros();

  for (int r = 0; r < rows; ++r) {
    for (int c = 0; c < cols; ++c) {
      precision_t h_pos = input_electrodynamics[0].mlts_indices(r, c);
      precision_t v_pos = input_electrodynamics[0].lats_indices(r, c);
      int r_start, c_start;

      if (h_pos < 0)
        c_start = 0;
      else if (h_pos >= matToInterpolateOn.n_cols - 1)
        c_start = matToInterpolateOn.n_cols - 2;
      else
        c_start = h_pos;

      if (v_pos < 0)
        r_start = 0;
      else if (v_pos >= matToInterpolateOn.n_rows - 1)
        r_start = matToInterpolateOn.n_rows - 2;
      else
        r_start = v_pos;

      precision_t first_row_slope =
        matToInterpolateOn(r_start, c_start + 1) -
        matToInterpolateOn(r_start, c_start);
      precision_t second_row_slope =
        matToInterpolateOn(r_start + 1, c_start + 1) -
        matToInterpolateOn(r_start + 1, c_start);
      precision_t first_row_val =
        matToInterpolateOn(r_start, c_start) +
        first_row_slope * (h_pos - c_start);
      precision_t second_row_val =
        matToInterpolateOn(r_start + 1, c_start) +
        second_row_slope * (h_pos - c_start);
      //time to vertically linear interpolate these two values
      precision_t vertical_slope = second_row_val - first_row_val;
      precision_t final_value =
        first_row_val + vertical_slope * (v_pos - r_start);
      slice(r, c) = final_value;
    }
  }

  return slice;
}

// -----------------------------------------------------------------------------
// Get potential, electron energy flux, and electron average energy patterns
// -----------------------------------------------------------------------------

std::tuple<arma_cube,
    arma_mat,
    arma_mat> Electrodynamics::get_electrodynamics(arma_cube magLat,
arma_cube magLocalTime) {
  arma_cube pot;
  arma_mat eflux;
  arma_mat avee;

  if (!input_electrodynamics.empty()) {
    pot = get_potential(magLat, magLocalTime);
    eflux = get_eflux(magLat, magLocalTime);
    avee = get_avee(magLat, magLocalTime);
  } else {
    pot.set_size(magLat.n_rows, magLat.n_cols, magLat.n_slices);
    pot.zeros();
    eflux.set_size(magLat.n_rows, magLat.n_cols);
    eflux.zeros();
    avee.set_size(magLat.n_rows, magLat.n_cols);
    avee.ones();
  }

  return std::make_tuple(pot, eflux, avee);
}

// -----------------------------------------------------------------------------
// Check if the requested time is within bounds of the electrodynamics times.
// -----------------------------------------------------------------------------

bool Electrodynamics::check_times(double inputStartTime, double inputEndTime) {

  if (HaveElectrodynamicsFile) {
    std::vector<double> e_times = input_electrodynamics[0].times;
    int iLow = 0;
    int iHigh = e_times.size() - 1;
    return !(inputStartTime > e_times[iHigh] || inputEndTime < e_times[iLow]);
  } else
    return true;
}

// -----------------------------------------------------------------------------
// Set the user-defined time for the electrodynamics, and find the
// interpolation indices for use when the user requests parameters.
// -----------------------------------------------------------------------------

void Electrodynamics::set_time(double time) {
  std::string function = "Electrodynamics::set_time";
  static int iFunction = -1;
  report.enter(function, iFunction);

  time_needed = time;

  if (!HaveElectrodynamicsFile) {
    time_index = -1;
    report.exit(function);
    return;
  }

  int64_t iLow, iMid, iHigh, N;
  double interpolation_index, x, dt;
  double intime = time;

  std::vector<double> times = input_electrodynamics[0].times;
  // Check to see if the time is below the bottom time in the vector:
  iLow = 0;
  iHigh = times.size() - 1;

  if (intime < times[iLow]) {
    interpolation_index = 0.0;
    report.print(0, "Warning: current time below first available\n");
    report.print(0, "potential-vector time, using first time\n");
  } else if (intime > times[iHigh]) {
    interpolation_index = iHigh;
    report.print(0, "Warning: current time above last available\n");
    report.print(0, "potential-vector time, using last time\n");
  }

  // At this point, we know that it is somewhere between the highest
  // and lowest values:
  else {
    iMid = (iHigh + iLow) / 2;

    while (iHigh - iLow > 1) {
      // Break if iMid <= time < iMid+1
      if (times[iMid] == intime)
        break;

      if (times[iMid] <= intime &&
          times[iMid + 1] > intime)
        break;

      // Upper Half:
      if (times[iMid] < intime) {
        iLow = iMid;
        iMid = (iHigh + iLow) / 2;
      } else {
        iHigh = iMid;
        iMid = (iHigh + iLow) / 2;
      }
    }

    // At this point, time should be between iMid and iMid+1:

    dt = (times[iMid + 1] - times[iMid]);
    x = (intime - times[iMid]) / dt;

    interpolation_index = iMid + x;
  }

  time_index = interpolation_index;
  report.exit(function);
}

// -----------------------------------------------------------------------------
// Take the user defined lats/mlts and set the interpolation indices
// so the code can easy get the requested electrodynamics parameters later
// -----------------------------------------------------------------------------

void Electrodynamics::set_grid(arma_mat lats, arma_mat mlts) {
  std::string function = "Electrodynamics::set_grid";
  static int iFunction = -1;
  report.enter(function, iFunction);

  lats_needed = abs(lats);
  mlts_needed = mlts;

  //uses first input_electrodynamics struct
  arma_vec lat_search = input_electrodynamics[0].mlats;
  arma_vec mlt_search = input_electrodynamics[0].mlts;

  input_electrodynamics[0].lats_indices = get_interpolation_indices(lats_needed,
                                          lat_search);
  input_electrodynamics[0].mlts_indices = get_interpolation_indices(mlts_needed,
                                          mlt_search);

  // If we have read in a file, we need to create interpolation indices here:

  // This is a bit more complicated, since we need to loop through
  // all of the points in lats and use the 1d interpolation scheme to find
  // the index, store, and move to next point; then repeat with mlts.

  report.exit(function);
}

// -----------------------------------------------------------------------------
// Use binomial search to find location in vector and then get linear
// interpolation coefficients
// -----------------------------------------------------------------------------

arma_mat Electrodynamics::get_interpolation_indices(arma_mat vals,
                                                    arma_vec search) {
  arma_mat res(vals.n_rows, vals.n_cols, fill::zeros);

  for (int i = 0; i < vals.n_rows; ++i) {
    for (int j = 0; j < vals.n_cols; ++j) {

      precision_t in = vals(i, j);
      int64_t iLow, iMid, iHigh, N;
      double interpolation_index, x, dx;

      iLow = 0;
      iHigh = search.n_rows - 1;

      // Check to see if the time is below or above the vector:

      if (in < search(iLow) || in > search(iHigh))
        interpolation_index = 0.0;

      else {

        // At this point, we know that it is somewhere between the highest
        // and lowest values:

        iMid = (iHigh + iLow) / 2;

        while (iHigh - iLow > 1) {
          // Break if iMid <= time < iMid+1
          if (search[iMid] == in)
            break;

          if (search[iMid] <= in &&
              search[iMid + 1] > in)
            break;

          // Upper Half:
          if (search[iMid] < in) {
            iLow = iMid;
            iMid = (iHigh + iLow) / 2;
          } else {
            iHigh = iMid;
            iMid = (iHigh + iLow) / 2;
          }
        }

        // At this point, time should be between iMid and iMid+1:

        dx = (search[iMid + 1] - search[iMid]);
        x = (in - search[iMid]) / dx;

        interpolation_index = iMid + x;
      }

      res(i, j) = interpolation_index;
    }
  }

  return res;
}

// -----------------------------------------------------------------------------
// Functions to set indices for models that need them
// -----------------------------------------------------------------------------

void Electrodynamics::set_imf_bx(precision_t value) {
  imf_bx_needed = value;
}

void Electrodynamics::set_imf_by(precision_t value) {
  imf_by_needed = value;
}

void Electrodynamics::set_imf_bz(precision_t value) {
  imf_bz_needed = value;
}

void Electrodynamics::set_sw_v(precision_t value) {
  sw_v_needed = value;
}

void Electrodynamics::set_sw_n(precision_t value) {
  sw_n_needed = value;
}

void Electrodynamics::set_hp(precision_t value) {
  hp_needed = value;
}

void Electrodynamics::set_au(precision_t value) {
  au_needed = value;
}

void Electrodynamics::set_al(precision_t value) {
  al_needed = value;
}

void Electrodynamics::set_ae(precision_t value) {
  ae_needed = value;
}

void Electrodynamics::set_kp(precision_t value) {
  kp_needed = value;
}

// --------------------------------------------------------------------------
// check to see if class is ok
// --------------------------------------------------------------------------

bool Electrodynamics::is_ok() {
  return IsOk;
}

/*
  Copyright (c) 2026 Atsuto Seko

  Permission is hereby granted, free of charge, to any person obtaining
  a copy of this software and associated documentation files (the
  "Software"), to deal in the Software without restriction, including
  without limitation the rights to use, copy, modify, merge, publish,
  distribute, sublicense, and/or sell copies of the Software, and to
  permit persons to whom the Software is furnished to do so, subject to
  the following conditions:

  The above copyright notice and this permission notice shall be
  included in all copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
  LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
  OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include <cstdlib>
#include <string>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>
#include <fstream>
#include <sstream>
#include <typeinfo>

#include "KIM_LogMacros.hpp"
#include "KIM_ModelDriverHeaders.hpp"

#include "polymlp/polymlp_mlpcpp.h"
#include "polymlp_kim.h"
#include "ndarray.hpp"

using namespace std;
using namespace model_driver_Tersoff;
// using namespace model_driver_polymlp;


extern "C" {
  // Can't be both static and extern.  But this is not needed: extern
  // just tells the compiler to use C naming conventions for the
  // function instead of C++ (name mangling).  All functions except
  // for the dirver creation are not linked at compile time but passed
  // at runtime by a pointer, so they do not need to be extern (and
  // should not be, there could be namespace clashes).
  int model_driver_create(KIM::ModelDriverCreate * const,
                          const KIM::LengthUnit,
                          const KIM::EnergyUnit,
                          const KIM::ChargeUnit,
                          const KIM::TemperatureUnit,
                          const KIM::TimeUnit);
}

// For some reason, KIM wants a pointer to that instead of just
// copying a bool. So be it, I'll define a constant, which is false
// since we do use ghost particles' neighbors.
static const int doesnt_use_ghost_neighbors = 0;

// LOCAL DEFINITIONS ///////////////////////////////////////////////////

// Helper to trim a string. For some reason C++ doesn't provide this.
/*
static string trim(const string &s)
{
    string::const_iterator it = s.begin();
    while (it != s.end() && isspace(*it))
        it++;

    string::const_reverse_iterator rit = s.rbegin();
    while (rit.base() != it && isspace(*rit))
        rit++;

    return string(it, rit.base());
}
*/


// WRAPPERS AND INTERFACE TO KIM ///////////////////////////////////////

static int 
compute_arguments_create(
    const KIM::ModelCompute * const, // unused
    KIM::ModelComputeArgumentsCreate * const model_compute_arguments_create) 
{
  int error =

    // Tell KIM what inputs/outputs are supported and how.
    model_compute_arguments_create->SetArgumentSupportStatus(
      KIM::COMPUTE_ARGUMENT_NAME::partialEnergy,
      KIM::SUPPORT_STATUS::optional)
    ||
    model_compute_arguments_create->SetArgumentSupportStatus(
      KIM::COMPUTE_ARGUMENT_NAME::partialParticleEnergy,
      KIM::SUPPORT_STATUS::optional)
    ||
    model_compute_arguments_create->SetArgumentSupportStatus(
      KIM::COMPUTE_ARGUMENT_NAME::partialForces,
      KIM::SUPPORT_STATUS::optional)
    ||
    model_compute_arguments_create->SetArgumentSupportStatus(
      KIM::COMPUTE_ARGUMENT_NAME::partialVirial,
      KIM::SUPPORT_STATUS::optional)
    ||
    model_compute_arguments_create->SetArgumentSupportStatus(
      KIM::COMPUTE_ARGUMENT_NAME::partialParticleVirial,
      KIM::SUPPORT_STATUS::optional)

    // Register callback support.
    ||
    model_compute_arguments_create->SetCallbackSupportStatus(
      KIM::COMPUTE_CALLBACK_NAME::ProcessDEDrTerm,
      KIM::SUPPORT_STATUS::optional)
    ||
    model_compute_arguments_create->SetCallbackSupportStatus(
      KIM::COMPUTE_CALLBACK_NAME::ProcessD2EDr2Term,
      KIM::SUPPORT_STATUS::notSupported)
    ;

  return error;
}
#undef KIM_LOGGER_OBJECT_NAME

#define KIM_LOGGER_OBJECT_NAME model_compute

static int 
compute(
    const KIM::ModelCompute * const model_compute,
    const KIM::ModelComputeArguments * const model_compute_arguments) 
{
    PolymlpKIM* polymlp_kim;
    model_compute->GetModelBufferPointer(reinterpret_cast<void **>(&polymlp_kim));

    // Unpack data.
    const int * n_atoms;
    const int * atom_types;
    const int * contributing;
    const double * atom_coords_ptr;

    double * energy;
    double * atom_energy;
    double * forces_ptr;
    double * virial;
    double * particle_virial_ptr;

    int error =

        // Input
        model_compute_arguments->GetArgumentPointer(
          KIM::COMPUTE_ARGUMENT_NAME::numberOfParticles, &n_atoms)
        ||
        model_compute_arguments->GetArgumentPointer(
          KIM::COMPUTE_ARGUMENT_NAME::particleSpeciesCodes, &atom_types)
        ||
        model_compute_arguments->GetArgumentPointer(
          KIM::COMPUTE_ARGUMENT_NAME::particleContributing, &contributing)
        ||
        model_compute_arguments->GetArgumentPointer(
          KIM::COMPUTE_ARGUMENT_NAME::coordinates, &atom_coords_ptr)

        // Output
        ||
        model_compute_arguments->GetArgumentPointer(
          KIM::COMPUTE_ARGUMENT_NAME::partialEnergy, &energy)
        ||
        model_compute_arguments->GetArgumentPointer(
          KIM::COMPUTE_ARGUMENT_NAME::partialParticleEnergy, &atom_energy)
        ||
        model_compute_arguments->GetArgumentPointer(
          KIM::COMPUTE_ARGUMENT_NAME::partialForces, &forces_ptr)
        ||
        model_compute_arguments->GetArgumentPointer(
          KIM::COMPUTE_ARGUMENT_NAME::partialVirial, &virial)
        ||
        model_compute_arguments->GetArgumentPointer(
          KIM::COMPUTE_ARGUMENT_NAME::partialParticleVirial, &particle_virial_ptr);

    if (error) return error;

    int compute_process_dEdr;
    error =
        model_compute_arguments->IsCallbackPresent(
            KIM::COMPUTE_CALLBACK_NAME::ProcessDEDrTerm, &compute_process_dEdr);

    if (error) return error;

    // Wrap some stuff for convenience.
    Array2D<const double> atom_coords(atom_coords_ptr, *n_atoms, 3);
    Array2D<double> f(forces_ptr, *n_atoms, 3);
    Array2D<double>* forces = forces_ptr ? &f : NULL;
    Array2D<double> v(particle_virial_ptr, *n_atoms, 6);
    Array2D<double>* particle_virial = particle_virial_ptr ? &v : NULL;

    // Do the compute.
    try {
      polymlp_kim.compute(
          *model_compute_arguments,
          *n_atoms,
          atom_types,
          contributing,
          atom_coords,
          energy,
          atom_energy,
          forces,
          virial,
          particle_virial,
          compute_process_dEdr);
    } catch (const exception& e) {
      LOG_ERROR(string("compute: ") + e.what());
      return 1;
    }

    return 0;
}

#undef KIM_LOGGER_OBJECT_NAME


static int 
compute_arguments_destroy(
    const KIM::ModelCompute * const, // all ununsed
    KIM::ModelComputeArgumentsDestroy * const) {
    // We did not allocate anything for compute_arguments_create(), thus
    // no cleanup is needed.

    return 0;
}


#define KIM_LOGGER_OBJECT_NAME model_destroy

static int destroy(KIM::ModelDestroy * const model_destroy) {
    PolymlpKIM* polymlp_kim;
    model_destroy->GetModelBufferPointer(reinterpret_cast<void **>(&polymlp_kim));

    if (polymlp_kim != NULL) {
        delete polymlp_kim;
    } 
    else {
        LOG_ERROR("destroy: tried to destroy a model driver that is already null");
    }
    return 0;
}
#undef KIM_LOGGER_OBJECT_NAME


#define KIM_LOGGER_OBJECT_NAME model_driver_create

static int
init_unit_conv(
    KIM::ModelDriverCreate * const model_driver_create,
    const KIM::LengthUnit length_unit,
    const KIM::EnergyUnit energy_unit,
    const KIM::ChargeUnit charge_unit,
    const KIM::TemperatureUnit temperature_unit,
    const KIM::TimeUnit time_unit,
    double& length_conv,
    double& inv_length_conv,
    double& energy_conv,
    double& inv_energy_conv,
    double& charge_conv) {
  int error;

  // Length ////////////////////////////////////////////////////////////
  error = model_driver_create->ConvertUnit(KIM::LENGTH_UNIT::A,
                                           KIM::ENERGY_UNIT::eV,
                                           KIM::CHARGE_UNIT::e,
                                           KIM::TEMPERATURE_UNIT::K,
                                           KIM::TIME_UNIT::ps,
                                           length_unit, energy_unit,
                                           charge_unit,
                                           temperature_unit, time_unit,
                                           1.0, 0.0, 0.0, 0.0, 0.0,
                                           &length_conv);
  if (error) {
    LOG_ERROR("Error returned by KIM's ConvertUnit() when trying to "
              "get length units.");
    return error;
  }

  // Inverse length ////////////////////////////////////////////////////
  error = model_driver_create->ConvertUnit(KIM::LENGTH_UNIT::A,
                                           KIM::ENERGY_UNIT::eV,
                                           KIM::CHARGE_UNIT::e,
                                           KIM::TEMPERATURE_UNIT::K,
                                           KIM::TIME_UNIT::ps,
                                           length_unit, energy_unit,
                                           charge_unit,
                                           temperature_unit, time_unit,
                                           -1.0, 0.0, 0.0, 0.0, 0.0,
                                           &inv_length_conv);
  if (error) {
    LOG_ERROR("Error returned by KIM's ConvertUnit() when trying to "
              "get inverse length units.");
    return error;
  }

  // Energy ////////////////////////////////////////////////////////////
  error = model_driver_create->ConvertUnit(KIM::LENGTH_UNIT::A,
                                           KIM::ENERGY_UNIT::eV,
                                           KIM::CHARGE_UNIT::e,
                                           KIM::TEMPERATURE_UNIT::K,
                                           KIM::TIME_UNIT::ps,
                                           length_unit, energy_unit,
                                           charge_unit,
                                           temperature_unit, time_unit,
                                           0.0, 1.0, 0.0, 0.0, 0.0,
                                           &energy_conv);
  if (error) {
    LOG_ERROR("Error returned by KIM's ConvertUnit() when trying to "
              "get energy units.");
    return error;
  }

  // Inverse energy ////////////////////////////////////////////////////
  error = model_driver_create->ConvertUnit(KIM::LENGTH_UNIT::A,
                                           KIM::ENERGY_UNIT::eV,
                                           KIM::CHARGE_UNIT::e,
                                           KIM::TEMPERATURE_UNIT::K,
                                           KIM::TIME_UNIT::ps,
                                           length_unit, energy_unit,
                                           charge_unit,
                                           temperature_unit, time_unit,
                                           0.0, -1.0, 0.0, 0.0, 0.0,
                                           &inv_energy_conv);
  if (error) {
    LOG_ERROR("Error returned by KIM's ConvertUnit() when trying to "
              "get inverse energy units.");
    return error;
  }

  // Charge ////////////////////////////////////////////////////////////
  error = model_driver_create->ConvertUnit(KIM::LENGTH_UNIT::A,
                                           KIM::ENERGY_UNIT::eV,
                                           KIM::CHARGE_UNIT::e,
                                           KIM::TEMPERATURE_UNIT::K,
                                           KIM::TIME_UNIT::ps,
                                           length_unit, energy_unit,
                                           charge_unit,
                                           temperature_unit, time_unit,
                                           0.0, 0.0, 1.0, 0.0, 0.0,
                                           &charge_conv);
  if (error) {
    LOG_ERROR("Error returned by KIM's ConvertUnit() when trying to "
              "get charge units.");
    return error;
  }

  // KIM wants to know what happened, I guess. /////////////////////////
  error = model_driver_create->SetUnits(length_unit, energy_unit, charge_unit,
                                        KIM::TEMPERATURE_UNIT::unused,
                                        KIM::TIME_UNIT::unused);
  if (error) {
    LOG_ERROR("Error returned by KIM's SetUnits().");
    return error;
  }

  return 0;
}
#undef KIM_LOGGER_OBJECT_NAME


#define KIM_LOGGER_OBJECT_NAME model_driver_create

// For readability, finish_create() should follow
// model_driver_create(), so we define its interface here already.
static int
finish_create(KIM::ModelDriverCreate * const,
              const KIM::LengthUnit,
              const KIM::EnergyUnit,
              const KIM::ChargeUnit,
              const KIM::TemperatureUnit,
              const KIM::TimeUnit,
              const std::string&,
              const int,
              std::map<std::string, int>&);


int
model_driver_create(KIM::ModelDriverCreate * const model_driver_create,
                    const KIM::LengthUnit length_unit,
                    const KIM::EnergyUnit energy_unit,
                    const KIM::ChargeUnit charge_unit,
                    const KIM::TemperatureUnit temperature_unit,
                    const KIM::TimeUnit time_unit) {
  int error;

  // Get parameter files. //////////////////////////////////////////////
  int n_param_files;
  model_driver_create->GetNumberOfParameterFiles(&n_param_files);
  if (n_param_files != 1) {
    LOG_ERROR("This model driver requires exactly one parameter file")
    return 1;
  }

  const string * param_filename;
  error = model_driver_create->GetParameterFileName(0, &param_filename);
  if (error) {
    LOG_ERROR("Error returned by KIM's GetParameterFileName() "
              "for the first parameter file.");
    return 1;
  }

  // Get number and name of species. ///////////////////////////////////
  /*
  int n_spec = 0;
  std::map<std::string, int> type_map;
  PotentialVariant potential_variant;
  error =
    read_settings(model_driver_create, *settings_filename,
                  n_spec, type_map, potential_variant);
  if (error) {
    return error; // already logged.
  }
  */

  return finish_create(model_driver_create,
                       length_unit,
                       energy_unit,
                       charge_unit,
                       temperature_unit,
                       time_unit,
                       *param_filename,
                       n_spec,
                       type_map);
}

static int
finish_create(KIM::ModelDriverCreate * const model_driver_create,
              const KIM::LengthUnit length_unit,
              const KIM::EnergyUnit energy_unit,
              const KIM::ChargeUnit charge_unit,
              const KIM::TemperatureUnit temperature_unit,
              const KIM::TimeUnit time_unit,
              const string& param_filename,
              const int n_spec,
              map<string,int>& type_map) {
  int error;

  // Init unit conversion factors. /////////////////////////////////////
  // We are using LAMMPS "metal" units, i.e., eV and Ångstroms.
  double length_conv;
  double inv_length_conv;
  double energy_conv;
  double inv_energy_conv;
  double charge_conv;
  error = init_unit_conv(model_driver_create,
                         length_unit, 
                         energy_unit, 
                         charge_unit,
                         temperature_unit, 
                         time_unit,
                         length_conv,
                         inv_length_conv,
                         energy_conv,
                         inv_energy_conv,
                         charge_conv);
  if (error) {
    return error; // already logged.
  }


  // Init the core class. //////////////////////////////////////////////
  PolymlpKIM* polymlp_kim;
  try {
    polymlp_kim = new PolymlpKIM(
        param_filename, 
        energy_conv, 
        inv_energy_conv,
        length_conv, 
        inv_length_conv,
        charge_conv);
  } catch (const exception& e) {
    LOG_ERROR(string("model_driver_create: ") + e.what());
    return 1; // error
  }

  // Pass stuff to KIM. ////////////////////////////////////////////////
  model_driver_create->SetModelBufferPointer(static_cast<void *>(polymlp_kim));
  // TODO: How to set cutoff_ptr.
  model_driver_create->SetInfluenceDistancePointer(tersoff->cutoff_ptr());
  model_driver_create->SetNeighborListPointers(1, tersoff->cutoff_ptr(),
                                               &doesnt_use_ghost_neighbors);
  error = model_driver_create->SetModelNumbering(KIM::NUMBERING::zeroBased);
  if (error) {
    LOG_ERROR("Error returned by KIM's SetModelNumbering().");
    delete tersoff;
    return 1;
  }

  // Register parameters.
  /*
  error = reg_params(model_driver_create, tersoff);
  if (error) {
    delete tersoff;
    return error; // logging already done in reg_params()
  }
  */

  // Use function pointer definitions to statically verify correct prototypes.
  KIM::ModelComputeArgumentsCreateFunction * kim_ca_create
    = &compute_arguments_create;
  KIM::ModelComputeFunction * kim_compute = &compute<T>;
  // KIM::ModelRefreshFunction * kim_refresh = &refresh<T>;
  // KIM::ModelWriteParameterizedModelFunction * kim_write_params
  //   = &write_parameterized_model<T>;
  KIM::ModelComputeArgumentsDestroyFunction * kim_ca_destroy
    = &compute_arguments_destroy;
  KIM::ModelDestroyFunction * kim_destroy = &destroy;

  // Register the function pointers.
  error =
    model_driver_create->SetRoutinePointer(
      KIM::MODEL_ROUTINE_NAME::ComputeArgumentsCreate,
      KIM::LANGUAGE_NAME::cpp, true,
      reinterpret_cast<KIM::Function *>(kim_ca_create))
    ||
    model_driver_create->SetRoutinePointer(
      KIM::MODEL_ROUTINE_NAME::Compute,
      KIM::LANGUAGE_NAME::cpp, true,
      reinterpret_cast<KIM::Function *>(kim_compute))
    ||
    /*
    model_driver_create->SetRoutinePointer(
      KIM::MODEL_ROUTINE_NAME::Refresh,
      KIM::LANGUAGE_NAME::cpp, false,
      reinterpret_cast<KIM::Function *>(kim_refresh))
    ||
    model_driver_create->SetRoutinePointer(
      KIM::MODEL_ROUTINE_NAME::WriteParameterizedModel,
      KIM::LANGUAGE_NAME::cpp, false,
      reinterpret_cast<KIM::Function *>(kim_write_params))
    ||
    */
    model_driver_create->SetRoutinePointer(
      KIM::MODEL_ROUTINE_NAME::ComputeArgumentsDestroy,
      KIM::LANGUAGE_NAME::cpp, true,
      reinterpret_cast<KIM::Function *>(kim_ca_destroy))
    ||
    model_driver_create->SetRoutinePointer(
      KIM::MODEL_ROUTINE_NAME::Destroy,
      KIM::LANGUAGE_NAME::cpp, true,
      reinterpret_cast<KIM::Function *>(kim_destroy));
  if (error) {
    LOG_ERROR("Error returned by KIM's SetRoutinePointer().");
    delete tersoff;
    return 1;
  }

  return 0;
}

#undef KIM_LOGGER_OBJECT_NAME

// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All right reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Radu Serban
// =============================================================================
//
// M113 suspension subsystem.
//
// =============================================================================

#include "m113/M113_Suspension.h"

using namespace chrono;
using namespace chrono::vehicle;

namespace m113 {

// -----------------------------------------------------------------------------
// Static variables
// -----------------------------------------------------------------------------
const double M113_Suspension::m_arm_mass = 1;
const ChVector<> M113_Suspension::m_arm_inertia(1, 1, 1);
const double M113_Suspension::m_arm_radius = 0.2;

// -----------------------------------------------------------------------------
// M113 shock functor class - implements a (non)linear damper
// -----------------------------------------------------------------------------
class M113_ShockForce : public ChSpringForceCallback {
  public:
    M113_ShockForce() {
        //// TODO
    }

    virtual double operator()(double time, double rest_length, double length, double vel) {
        //// TODO
        return 0;
    }
};

// -----------------------------------------------------------------------------
// M113 torsion-bar force - implements a (non)linear rotational elastic force
// -----------------------------------------------------------------------------
class M113_TorsionForce : public ChTorsionForce {
  public:
    M113_TorsionForce() {
        //// TODO
    }
};

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
M113_Suspension::M113_Suspension() : ChLinearDamperRWAssembly("M113_Suspension") {
    m_shock_forceCB = new M113_ShockForce();
    m_torsion_force = new M113_TorsionForce();
}

M113_Suspension::~M113_Suspension() {
    delete m_shock_forceCB;
    delete m_torsion_force;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
const ChVector<> M113_Suspension::GetLocation(PointId which) {
    switch (which) {
        case ARM:
            return ChVector<>(0, 0, 0);
        case ARM_CHASSIS:
            return ChVector<>(0, 0, 0);
        case SHOCK_A:
            return ChVector<>(0, 0, 0);
        case SHOCK_C:
            return ChVector<>(0, 0, 0);
        default:
            return ChVector<>(0, 0, 0);
    }
}

}  // end namespace m113
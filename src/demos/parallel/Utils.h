// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Milad Rakhsha
// =============================================================================

#include <cstdio>
#include <vector>
#include <cmath>
#include <memory>
#include "chrono_parallel/physics/ChSystemParallel.h"
#include "chrono/ChConfig.h"
#include "chrono/utils/ChUtilsCreators.h"
#include "chrono/utils/ChUtilsInputOutput.h"
#include "chrono/utils/ChUtilsGeometry.h"
#include "chrono/utils/ChUtilsGenerators.h"
#include "chrono_parallel/physics/Ch3DOFContainer.h"

using namespace chrono;
using namespace chrono::collision;

void OutputBlazeMatrix(CompressedMatrix<real> src, std::string filename) {
    const char* numformat = "%.16g";
    ChStreamOutAsciiFile stream(filename.c_str());
    stream.SetNumFormat(numformat);
    for (int i = 0; i < src.rows(); ++i) {
        for (CompressedMatrix<real>::Iterator it = src.begin(i); it != src.end(i); ++it) {
            stream << it->value() << ",";
        }
        stream << "\n";
    }
}

void OutputBlazeVector(DynamicVector<real> src, std::string filename) {
    const char* numformat = "%.16g";
    ChStreamOutAsciiFile stream(filename.c_str());
    stream.SetNumFormat(numformat);

    for (int i = 0; i < src.size(); i++)
        stream << src[i] << ",";
}

void writeCSV(ChSystemParallel* msystem, int out_frame, std::string out_folder) {
    char filename2[100];
    sprintf(filename2, "%s/data_%d.csv", out_folder.c_str(), out_frame + 1);

    const std::string& delim = ",";
    utils::CSV_writer csv(delim);
    int numMarkers = msystem->data_manager->host_data.pos_rigid.size();
    csv << "t,x,y,z,vx,vy,vz,|U|" << std::endl;
    for (int i = 0; i < numMarkers; i++) {
        real3 pos3 = msystem->data_manager->host_data.pos_rigid[i];
        real vx = msystem->data_manager->host_data.v[6 * i];
        real vy = msystem->data_manager->host_data.v[6 * i + 1];
        real vz = msystem->data_manager->host_data.v[6 * i + 2];
        real u = sqrt(vx * vx + vy * vy + vz * vz);
        csv << msystem->GetChTime() << pos3.x << pos3.y << pos3.z << vx << vy << vz << std::endl;
    }

    csv.write_to_file(filename2);
}

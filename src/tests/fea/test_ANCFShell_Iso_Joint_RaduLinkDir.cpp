#include "chrono/physics/ChSystem.h"
#include "chrono/lcp/ChLcpIterativeMINRES.h"
#include "chrono_fea/ChElementShellANCF.h"
#include "chrono_fea/ChMesh.h"
#include "chrono_fea/ChLinkPointFrame.h"
#include "chrono_fea/ChLinkDirFrame.h"
#include "chrono/utils/ChUtilsInputOutput.h"
#include "chrono/utils/ChUtilsValidation.h"
#include "chrono_mkl/ChLcpMklSolver.h"
#include "chrono/physics/ChBody.h"

////#include <float.h>
////unsigned int fp_control_state = _controlfp(_EM_INEXACT, _MCW_EM);

using namespace chrono;
using namespace fea;

// ========================================================================

double time_step = 0.001;  // Time step
int num_steps = 5000;      // Number of time steps for unit test (range 1 to 4000)


// Bodies

bool include_bodies = true;
ChSharedPtr<ChBody> BGround;
ChSharedPtr<ChBody> Body_2;
ChSharedPtr<ChBody> Body_3;

// Mesh
bool include_mesh = true;
ChSharedPtr<ChMesh> my_mesh;
ChSharedPtr<ChNodeFEAxyzD> NodeFirst;
ChSharedPtr<ChNodeFEAxyzD> Node2;
ChSharedPtr<ChNodeFEAxyzD> Node3;
ChSharedPtr<ChNodeFEAxyzD> Node4;

// Joints
bool include_joints = true;
ChSharedPtr<ChLinkLockLock> my_link_12;
ChSharedPtr<ChLinkLockRevolute> my_link_23;

// Body-mesh constraints
bool include_constraints = true;
ChSharedPtr<ChLinkPointFrame> constraint_hinge;
ChSharedPtr<ChLinkDirFrame> constraintDir;

// Output data
bool print_data = true;
bool save_data = false;
utils::Data m_data;

// ========================================================================

void AddBodies(ChSystem& my_system) {
    if (!include_bodies)
        return;

    // Defining the Body 1
    BGround = ChSharedPtr<ChBody>(new ChBody);
    my_system.AddBody(BGround);
    BGround->SetIdentifier(1);
    BGround->SetBodyFixed(true);
    BGround->SetCollide(false);
    BGround->SetMass(1);
    BGround->SetInertiaXX(ChVector<>(1, 1, 0.2));
    BGround->SetPos(ChVector<>(-2, 0, 0));  // Y = -1m
    ChQuaternion<> rot = Q_from_AngX(0.0);
    BGround->SetRot(rot);
    // Defining the Body 2
    Body_2 = ChSharedPtr<ChBody>(new ChBody);
    my_system.AddBody(Body_2);
    Body_2->SetIdentifier(2);
    Body_2->SetBodyFixed(false);
    Body_2->SetCollide(false);
    Body_2->SetMass(1);
    Body_2->SetInertiaXX(ChVector<>(1, 1, 0.2));
    Body_2->SetPos(ChVector<>(-1, 0, 0));  // Y = -1m
    Body_2->SetRot(rot);
    // Defining the Body 3
    Body_3 = ChSharedPtr<ChBody>(new ChBody);
    my_system.AddBody(Body_3);
    Body_3->SetIdentifier(3);
    Body_3->SetBodyFixed(false);
    Body_3->SetCollide(false);
    Body_3->SetMass(2);
    Body_3->SetInertiaXX(ChVector<>(2, 2, 0.4));
    Body_3->SetPos(ChVector<>(0.25, 0, 0));
    Body_3->SetRot(rot);
}

// ========================================================================

void AddMesh(ChSystem& my_system) {
    if (!include_mesh)
        return;

    // Create a mesh, that is a container for groups of elements and their referenced nodes.
    my_mesh = ChSharedPtr<ChMesh>(new ChMesh);
    // Geometry of the plate
    double plate_lenght_x = 1;
    double plate_lenght_y = 1;
    double plate_lenght_z = 0.01;  // small thickness
    // Specification of the mesh
    const int numDiv_x = 1;
    const int numDiv_y = 1;
    const int numDiv_z = 1;
    const int N_x = numDiv_x + 1;
    const int N_y = numDiv_y + 1;
    const int N_z = numDiv_z + 1;
    // Number of elements in the z direction is considered as 1
    int TotalNumElements = numDiv_x * numDiv_y;
    int TotalNumNodes = (numDiv_x + 1) * (numDiv_y + 1);
    // For uniform mesh
    double dx = plate_lenght_x / numDiv_x;
    double dy = plate_lenght_y / numDiv_y;
    double dz = plate_lenght_z / numDiv_z;
    int MaxMNUM = 0;
    int MTYPE = 0;
    int MaxLayNum = 0;
    ChMatrixDynamic<double> COORDFlex(TotalNumNodes, 6);
    ChMatrixDynamic<double> VELCYFlex(TotalNumNodes, 6);
    ChMatrixDynamic<int> NumNodes(TotalNumElements, 4);
    ChMatrixDynamic<int> LayNum(TotalNumElements, 1);  // Only one layer in this unit test
    ChMatrixDynamic<int> NDR(TotalNumNodes, 6);
    ChMatrixDynamic<double> ElemLengthXY(TotalNumElements, 2);
    ChMatrixNM<double, 10, 12> MPROP;
    ChMatrixNM<int, 10, 7> MNUM;
    ChMatrixNM<int, 10, 1> NumLayer;
    double LayPROP[10][7][2];

    //!------------------------------------------------!
    //!--------------- Element data--------------------!
    //!------------------------------------------------!
    for (int i = 0; i < TotalNumElements; i++) {
        // All the elements belong to the same layer, e.g layer number 1.
        LayNum(i, 0) = 1;
        // Node number of the 4 nodes which creates element i.
        // The nodes are distributed this way. First in the x direction for constant
        // y when max x is reached go to the next level for y by doing the same
        // distribution but for y+1 and keep doing until y limit is reached. Node
        // number start from 1.
        NumNodes(i, 0) = (i / (numDiv_x)) * (N_x) + i % numDiv_x;
        NumNodes(i, 1) = (i / (numDiv_x)) * (N_x) + i % numDiv_x + 1;
        NumNodes(i, 2) = (i / (numDiv_x)) * (N_x) + i % numDiv_x + N_x;
        NumNodes(i, 3) = (i / (numDiv_x)) * (N_x) + i % numDiv_x + 1 + N_x;
        // Let's keep the element length a fixed number in both direction. (uniform
        // distribution of nodes in both direction)
        // If the user interested in non-uniform mesh it can be manipulated here.
        ElemLengthXY(i, 0) = dx;
        ElemLengthXY(i, 1) = dy;
        if (MaxLayNum < LayNum(i, 0)) {
            MaxLayNum = LayNum(i, 0);
        }
    }

    //!----------------------------------------------!
    //!--------- NDR,COORDFlex,VELCYFlex-------------!
    //!----------------------------------------------!
    for (int i = 0; i < TotalNumNodes; i++) {
        // If the node is the first node from the left side fix the x,y,z degree of
        // freedom+ d/dx, d/dy,d/dz as well. 1for constrained 0 for ...
        //-The NDR array is used to define the degree of freedoms that are
        // constrained in the 6 variable explained above.
        NDR(i, 0) = (i % (numDiv_x + 1) == 0) ? 1 : 0;
        NDR(i, 1) = (i % (numDiv_x + 1) == 0) ? 1 : 0;
        NDR(i, 2) = (i % (numDiv_x + 1) == 0) ? 1 : 0;
        NDR(i, 3) = (i % (numDiv_x + 1) == 0) ? 1 : 0;
        NDR(i, 4) = (i % (numDiv_x + 1) == 0) ? 1 : 0;
        NDR(i, 5) = (i % (numDiv_x + 1) == 0) ? 1 : 0;

        //-COORDFlex are the initial coordinates for each node,
        // the first three are the position and the last three are the slope
        // coordinates for the tangent vector.
        // Note that i starts from 0 but nodes starts from 1
        COORDFlex(i, 0) = (i % (numDiv_x + 1)) * dx + 0.5;
        COORDFlex(i, 1) = (i / (numDiv_x + 1)) % (numDiv_y + 1) * dy;
        COORDFlex(i, 2) = (i) / ((numDiv_x + 1) * (numDiv_y + 1)) * dz;
        COORDFlex(i, 3) = 0;
        COORDFlex(i, 4) = 0;
        COORDFlex(i, 5) = 1;
        //-VELCYFlex is essentially the same as COORDFlex, but for the initial
        // velocity instead of position.
        // let's assume zero initial velocity for nodes
        VELCYFlex(i, 0) = 0;
        VELCYFlex(i, 1) = 0;
        VELCYFlex(i, 2) = 0;
        VELCYFlex(i, 3) = 0;
        VELCYFlex(i, 4) = 0;
        VELCYFlex(i, 5) = 0;
    }

    //!------------------------------------------------!
    //!------------- Read Layer Data-------------------!
    //!------------------------------------------------!
    for (int i = 0; i < MaxLayNum; i++) {
        NumLayer(i, 0) = i + 1;
        // For each layer define some properties
        // For multilayer problem the user should modify the following loop as they
        // would like
        for (int j = 0; j < NumLayer(i, 0); j++) {
            LayPROP[i][j][0] = dz;  // Layerheight
            LayPROP[i][j][1] = 0;   // PlyAngle
            MNUM[i][j] = 1;         // Material_ID
            if (MaxMNUM < MNUM(i, j))
                MaxMNUM = MNUM(i, j);
        }
    }

    //!------------------------------------------------!
    //!------------ Input Material Data----------------!
    //!------------------------------------------------!
    for (int i = 0; i < MaxMNUM; i++) {
        double nu_coef = 0.3;
        MTYPE = 2;  // The user must use orthotropic input (MTYPE=2) to introduce isotropic material
        // properties for this unit test

        if (MTYPE == 2) {
            MPROP(i, 0) = 500;      // Density [kg/m3]
            MPROP(i, 1) = 2.1E+08;  // Ex//
            MPROP(i, 2) = 2.1E+08;  // Ey
            MPROP(i, 3) = 2.1E+08;  // Ez
            // Additional information for the Type 2 of Material.
            MPROP(i, 4) = 0.3;                                    // nuxy
            MPROP(i, 5) = 0.3;                                    // nuxz
            MPROP(i, 6) = 0.3;                                    // nuyz
            MPROP(i, 7) = MPROP(i, 1) / 2.0 / (1 + MPROP(i, 6));  // Gxy
            MPROP(i, 8) = MPROP(i, 1) / 2.0 / (1 + MPROP(i, 6));  // Gxz
            MPROP(i, 9) = MPROP(i, 1) / 2.0 / (1 + MPROP(i, 6));  // Gyz
        }
    }
    ChSharedPtr<ChContinuumElastic> mmaterial(new ChContinuumElastic);
    // Adding the nodes to the mesh
    int i = 0;

    while (i < TotalNumNodes) {
        ChSharedPtr<ChNodeFEAxyzD> node(
            new ChNodeFEAxyzD(ChVector<>(COORDFlex(i, 0), COORDFlex(i, 1), COORDFlex(i, 2)),
                              ChVector<>(COORDFlex(i, 3), COORDFlex(i, 4), COORDFlex(i, 5))));
        node->SetMass(0.0);
        my_mesh->AddNode(node);
        if (NDR(i, 0) == 1 && NDR(i, 1) == 1 && NDR(i, 2) == 1 && NDR(i, 3) == 1 && NDR(i, 4) == 1 && NDR(i, 5) == 1) {
            // node->SetFixed(true);
        }
        i++;
    }

    NodeFirst = ChSharedPtr<ChNodeFEAxyzD>(my_mesh->GetNode(0).DynamicCastTo<ChNodeFEAxyzD>());
    Node2 = ChSharedPtr<ChNodeFEAxyzD>(my_mesh->GetNode(1).DynamicCastTo<ChNodeFEAxyzD>());
    Node3 = ChSharedPtr<ChNodeFEAxyzD>(my_mesh->GetNode(2).DynamicCastTo<ChNodeFEAxyzD>());
    Node4 = ChSharedPtr<ChNodeFEAxyzD>(my_mesh->GetNode(3).DynamicCastTo<ChNodeFEAxyzD>());

    int elemcount = 0;
    while (elemcount < TotalNumElements) {
        ChSharedPtr<ChElementShellANCF> element(new ChElementShellANCF);
        // Save material data into InertFlexVec(98x1) at each layer
        ChMatrixNM<double, 98, 1> InertFlexVec;
        InertFlexVec.Reset();
        double TotalThickness;  // element thickness
        TotalThickness = 0.0;
        int i = elemcount;
        for (int j = 0; j < NumLayer(LayNum(i, 0) - 1, 0); j++) {
            int ij = 14 * j;
            InertFlexVec(ij) = MPROP[MNUM[LayNum(i, 0) - 1][j] - 1][0];  // Density
            InertFlexVec(ij + 1) = ElemLengthXY(i, 0);                   // EL
            InertFlexVec(ij + 2) = ElemLengthXY(i, 1);                   // EW
            InertFlexVec(ij + 3) = LayPROP[LayNum(i, 0) - 1][j][0];      // Thickness per layer
            TotalThickness += InertFlexVec(ij + 3);
            InertFlexVec(ij + 4) = LayPROP[LayNum(i, 0) - 1][j][1];           // Fiber angle
            InertFlexVec(ij + 5) = MPROP[MNUM[LayNum(i, 0) - 1][j] - 1][1];   // Ex
            InertFlexVec(ij + 6) = MPROP[MNUM[LayNum(i, 0) - 1][j] - 1][2];   // Ey
            InertFlexVec(ij + 7) = MPROP[MNUM[LayNum(i, 0) - 1][j] - 1][3];   // Ez
            InertFlexVec(ij + 8) = MPROP[MNUM[LayNum(i, 0) - 1][j] - 1][4];   // nuxy
            InertFlexVec(ij + 9) = MPROP[MNUM[LayNum(i, 0) - 1][j] - 1][5];   // nuxz
            InertFlexVec(ij + 10) = MPROP[MNUM[LayNum(i, 0) - 1][j] - 1][6];  // nuyz
            InertFlexVec(ij + 11) = MPROP[MNUM[LayNum(i, 0) - 1][j] - 1][7];  // Gxy
            InertFlexVec(ij + 12) = MPROP[MNUM[LayNum(i, 0) - 1][j] - 1][8];  // Gxz
            InertFlexVec(ij + 13) = MPROP[MNUM[LayNum(i, 0) - 1][j] - 1][9];  // Gyz
        }
        ChMatrixNM<double, 7, 2> GaussZRange;
        GaussZRange.Reset();
        double CurrentHeight = 0.0;
        for (int j = 0; j < NumLayer(LayNum(i, 0) - 1, 0); j++) {
            double AA = (CurrentHeight / TotalThickness - 0.5) * 2.0;
            CurrentHeight += LayPROP[LayNum(i, 0) - 1][j][0];
            double AAA = (CurrentHeight / TotalThickness - 0.5) * 2.0;
            GaussZRange(j, 0) = AA;
            GaussZRange(j, 1) = AAA;
        }
        element->SetInertFlexVec(InertFlexVec);
        element->SetGaussZRange(GaussZRange);
        element->SetNodes(my_mesh->GetNode(NumNodes[elemcount][0]).DynamicCastTo<ChNodeFEAxyzD>(),
                          my_mesh->GetNode(NumNodes[elemcount][1]).DynamicCastTo<ChNodeFEAxyzD>(),
                          my_mesh->GetNode(NumNodes[elemcount][2]).DynamicCastTo<ChNodeFEAxyzD>(),
                          my_mesh->GetNode(NumNodes[elemcount][3]).DynamicCastTo<ChNodeFEAxyzD>());
        element->SetMaterial(mmaterial);
        element->SetNumLayers(NumLayer(LayNum(i, 0) - 1, 0));
        element->SetThickness(TotalThickness);
        element->SetElemNum(elemcount);
        element->SetAlphaDamp(0.08);
        element->Setdt(time_step);                 // dt to calculate DampingCoefficient
        element->SetGravityOn(false);              // turn gravity on/off
        ChMatrixNM<double, 35, 1> StockAlpha_EAS;  // StockAlpha(5*7,1): Max #Layer is 7
        StockAlpha_EAS.Reset();
        element->SetStockAlpha(StockAlpha_EAS);
        my_mesh->AddElement(element);
        elemcount++;
    }

    // Switch off mesh class gravity
    my_mesh->SetAutomaticGravity(false);
    // This is mandatory
    my_mesh->SetupInitial();
    // Remember to add the mesh to the system!
    my_system.Add(my_mesh);
}

// ========================================================================

void AddConstraints(ChSystem& my_system) {
    if (include_bodies && include_joints) {
        // Weld body_2 to ground body.
        my_link_12 = ChSharedPtr<ChLinkLockLock>(new ChLinkLockLock);
        my_link_12->Initialize(BGround, Body_2, ChCoordsys<>(ChVector<>(-2.0, 0, 0)));
        my_system.AddLink(my_link_12);

        // Add another revolute joint
        my_link_23 = ChSharedPtr<ChLinkLockRevolute>(new ChLinkLockRevolute);
        my_link_23->Initialize(Body_2, Body_3, ChCoordsys<>(ChVector<>(0, 0, 0), Q_from_AngX(CH_C_PI / 2.0)));
        my_system.AddLink(my_link_23);
    }

    if (include_bodies && include_mesh && include_constraints) {
        // Constraining a node to the truss
        constraint_hinge = ChSharedPtr<ChLinkPointFrame>(new ChLinkPointFrame);
        constraint_hinge->Initialize(NodeFirst, Body_3);
        my_system.Add(constraint_hinge);

        // This contraint means that rz will always be perpendicular to
        // the direction along the length of the link
        constraintDir = ChSharedPtr<ChLinkDirFrame>(new ChLinkDirFrame);
        constraintDir->Initialize(NodeFirst, Body_3);
        // constraintDir->SetDirectionInBodyCoords(ChVector<double>(0, 0, 1));
        constraintDir->SetDirectionInAbsoluteCoords(ChVector<double>(0, 0, 1));
        my_system.Add(constraintDir);
        GetLog() << constraintDir->GetDirection();
    }
}

// ========================================================================

void SaveData(ChSystem& my_system, utils::CSV_writer& csv, int it) {
    if (!save_data || !include_bodies || !include_mesh)
        return;

    m_data[0][it] = my_system.GetChTime();
    m_data[1][it] = Body_2->coord.pos.x;
    m_data[2][it] = Body_2->coord.pos.y;
    m_data[3][it] = Body_2->coord.pos.z;
    m_data[4][it] = Body_3->coord.pos.x;
    m_data[5][it] = Body_3->coord.pos.y;
    m_data[6][it] = Body_3->coord.pos.z;
    m_data[7][it] = NodeFirst->pos.x;
    m_data[8][it] = NodeFirst->pos.y;
    m_data[9][it] = NodeFirst->pos.z;

    m_data[10][it] = Node2->pos.x;
    m_data[11][it] = Node2->pos.y;
    m_data[12][it] = Node2->pos.z;
    m_data[13][it] = Node4->pos.x;
    m_data[14][it] = Node4->pos.y;
    m_data[15][it] = Node4->pos.z;
    m_data[16][it] = NodeFirst->D.x;
    m_data[17][it] = NodeFirst->D.y;
    m_data[18][it] = NodeFirst->D.z;
    csv << m_data[0][it] << m_data[1][it] << m_data[2][it] << m_data[3][it] << m_data[4][it] << m_data[5][it]
        << m_data[6][it] << m_data[7][it] << m_data[8][it] << m_data[9][it] << m_data[10][it] << m_data[11][it]
        << m_data[12][it] << m_data[13][it] << m_data[14][it] << m_data[15][it] << m_data[16][it] << m_data[17][it]
        << m_data[18][it] << std::endl;
}

// ========================================================================

int main(int argc, char* argv[]) {
    // Consistency
    include_joints = include_joints && include_bodies;
    include_constraints = include_constraints && include_bodies && include_mesh;

    // Definition of the model
    ChSystem my_system;

    my_system.Set_G_acc(ChVector<>(0, 0, -9.81));

    AddMesh(my_system);
    AddBodies(my_system);
    AddConstraints(my_system);

    // Set up linear solver
    ChLcpMklSolver* mkl_solver_stab = new ChLcpMklSolver;
    ChLcpMklSolver* mkl_solver_speed = new ChLcpMklSolver;
    my_system.ChangeLcpSolverStab(mkl_solver_stab);
    my_system.ChangeLcpSolverSpeed(mkl_solver_speed);
    mkl_solver_stab->SetSparsityPatternLock(true);
    mkl_solver_speed->SetSparsityPatternLock(true);
    my_system.Update();

    // Set up integrator
    my_system.SetIntegrationType(ChSystem::INT_HHT);
    ChSharedPtr<ChTimestepperHHT> mystepper = my_system.GetTimestepper().DynamicCastTo<ChTimestepperHHT>();
    mystepper->SetAlpha(-0.2);
    mystepper->SetMaxiters(1000);
    mystepper->SetTolerance(1e-07);
    mystepper->SetMode(ChTimestepperHHT::POSITION);
    mystepper->SetScaling(true);

    m_data.resize(19);
    for (size_t col = 0; col < 19; col++)
        m_data[col].resize(num_steps);
    utils::CSV_writer csv(" ");

    for (int it = 0; it < num_steps; it++) {
        my_system.DoStepDynamics(time_step);

        if (print_data) {
            std::cout << "\nTime t = " << my_system.GetChTime() << "s \n";

            if (include_bodies) {
                printf("Body_2 position: %12.4e  %12.4e  %12.4e\n", Body_2->coord.pos.x, Body_2->coord.pos.y,
                    Body_2->coord.pos.z);
                printf("Body_3 position: %12.4e  %12.4e  %12.4e\n", Body_3->coord.pos.x, Body_3->coord.pos.y,
                    Body_3->coord.pos.z);
                ChVector<> tip = Body_3->TransformPointLocalToParent(ChVector<>(0.25, 0, 0));
                printf("Body_3 tip:      %12.4e  %12.4e  %12.4e\n", tip.x, tip.y, tip.z);
            }

            if (include_mesh) {
                // std::cout << "nodetip->pos.z = " << Node4->pos.z << "\n";
                printf("Node position:   %12.4e  %12.4e  %12.4e\n", NodeFirst->pos.x, NodeFirst->pos.y, NodeFirst->pos.z);
                printf("Direction of node:  %12.4e  %12.4e  %12.4e\n", NodeFirst->D.x, NodeFirst->D.y, NodeFirst->D.z);
            }

            if (include_constraints) {
                // Get direction of constraint (in body local frame) and convert to global frame
                ChVector<> dirB = Body_3->TransformDirectionLocalToParent(constraintDir->GetDirection());
                printf("Direction on body:  %12.4e  %12.4e  %12.4e\n", dirB.x, dirB.y, dirB.z);
                // Direction along the body
                ChVector<> body_axis = Body_3->TransformDirectionLocalToParent(ChVector<>(0.25, 0, 0));
                printf("Body axis dir:      %12.4e  %12.4e  %12.4e\n", body_axis.x, body_axis.y, body_axis.z);
                // Body axis should always be perpendicular to node normal
                double dot = Vdot(body_axis, NodeFirst->D);
                printf("Dot product = %e\n", dot);

                ChMatrixNM<double, 3, 1> Cp = constraint_hinge->GetC();
                printf("Point constraint violations:      %12.4e  %12.4e  %12.4e\n", Cp.GetElement(0, 0), Cp.GetElement(1, 0), Cp.GetElement(2, 0));
                ChMatrixNM<double, 2, 1> Cd = constraintDir->GetC();
                printf("Direction constraint violations:  %12.4e  %12.4e\n", Cd.GetElement(0, 0), Cd.GetElement(1, 0));
            }

            if (include_joints) {
                ChMatrix<>* C12 = my_link_12->GetC();
                printf("Weld joint constraints: %12.4e  %12.4e  %12.4e  %12.4e  %12.4e  %12.4e\n", C12->GetElement(0, 0),
                    C12->GetElement(1, 0), C12->GetElement(2, 0), C12->GetElement(3, 0), C12->GetElement(4, 0),
                    C12->GetElement(5, 0));

                ChMatrix<>* C23 = my_link_23->GetC();
                printf("Rev joint constraints:  %12.4e  %12.4e  %12.4e  %12.4e  %12.4e\n", C23->GetElement(0, 0),
                    C23->GetElement(1, 0), C23->GetElement(2, 0), C23->GetElement(3, 0), C23->GetElement(4, 0));
            }
        }

        SaveData(my_system, csv, it);
    }

    csv.write_to_file("Body_Mesh_constraints.txt");

    return 0;
}

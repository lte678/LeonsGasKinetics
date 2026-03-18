// Example: Knudsen Pump
// Byambadorj et al. "Monolithic SOI through-wafer Knudsen pumps with mechanically robust Si channels". 2023
channelWidth = 1.45;
reservoirWidth = 10 + channelWidth;
reservoirHeight = reservoirWidth/2;
coldReservoirWidth = 6;
coldReservoirHeight = 3;

Point(1) = {-channelWidth/2, 0, 0, 1.0};
Point(2) = {channelWidth/2, 0, 0, 1.0};
Point(3) = {channelWidth/2, 12, 0, 1.0};
Point(4) = {-channelWidth/2, 12, 0, 1.0};
Point(9) = {-channelWidth/2, 14, 0, 1.0};
Point(10) = {channelWidth/2, 14, 0, 1.0};
Point(11) = {-reservoirWidth / 2.0, 14, 0, 1.0};
Point(12) = {reservoirWidth / 2.0, 14, 0, 1.0};
Point(13) = {reservoirWidth / 2.0, 14 + reservoirHeight, 0, 1.0};
Point(14) = {-reservoirWidth / 2.0, 14 + reservoirHeight, 0, 1.0};
Point(15) = {-coldReservoirWidth / 2.0, -coldReservoirHeight, 0, 1.0};
Point(16) = {coldReservoirWidth / 2.0, -coldReservoirHeight, 0, 1.0};
Point(17) = {coldReservoirWidth / 2.0, 0, 0, 1.0};
Point(18) = {-coldReservoirWidth / 2.0, 0, 0, 1.0};

Line(1) = {14, 13};
Line(2) = {13, 12};
Line(3) = {12, 10};
Line(4) = {10, 3};
Line(5) = {3, 2};
Line(6) = {2, 17};
Line(7) = {17, 16};
Line(8) = {16, 15};
Line(9) = {15, 18};
Line(10) = {18, 1};
Line(11) = {1, 4};
Line(12) = {4, 9};
Line(13) = {9, 11};
Line(14) = {11, 14};

Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14};

Plane Surface(1) = {-1};

MeshSize {11, 12, 13, 14} = 2.0
MeshSize {15, 16, 17, 18} = 1.0;
MeshSize {9, 4, 10, 3, 1, 2} = 0.5;

Mesh.Algorithm = 8;
Mesh 2;
Mesh.SubdivisionAlgorithm = 2;
RefineMesh;

Physical Curve("ChannelCold", 15) = {11, 5};
Physical Curve("ChannelHot", 16) = {12, 4};
Physical Curve("CavityHotL", 17) = {13};
Physical Curve("CavityHotR", 22) = {3};
Physical Curve("CavityCold", 18) = {10, 6};
Physical Curve("OpenCold", 19) = {8, 9, 7};
Physical Curve("OpenHot", 20) = {1, 14, 2};
Physical Surface("Symmetry", 21) = {1};

Mesh.MshFileVersion = 4.1;
Mesh.SaveAll = 1;
Save "knudsen_pump.msh";

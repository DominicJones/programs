Mesh.Format = 1;
Mesh.MshFileVersion = 2.2;

Mesh.MeshSizeMin = 1.0;
Mesh.MeshSizeMax = 1.0;
Mesh.MeshSizeFactor = 1.0;

lc = 0.2;
//+
SetFactory("OpenCASCADE");
Box(1) = {0, 0, 0, 1, 1, 1};
//+
Surface Loop(2) = {2, 3, 5, 1, 6, 4};
// Volume(30) = {2};
Physical Surface(31) = {2};
Physical Surface(32) = {3, 5, 1, 6, 4};
Physical Volume(100) = {1};

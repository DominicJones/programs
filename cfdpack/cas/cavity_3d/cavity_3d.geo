lc = 0.1;
Point(1) = {0, 0, 0, 1.0, lc};
Point(2) = {1, -0, 0, 1.0, lc};
Point(3) = {1, 1, 0, 1.0, lc};
Point(4) = {-0, 1, 0, 1.0, lc};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(5) = {3, 4, 1, 2};
Plane Surface(6) = {5};
Extrude {0, 0, 1} {
  Surface{6};
}
Surface Loop(29) = {28, 15, 6, 19, 23, 27};
Volume(30) = {29};
Physical Surface(31) = {15};
Physical Surface(32) = {19, 28, 6, 27, 23};
Physical Volume(100) = {30};

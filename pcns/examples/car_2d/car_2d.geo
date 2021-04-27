Mesh.Format = 1;
Mesh.MshFileVersion = 2.2;
//
// characteristic lengths [box, car]
//
box_cl = 1.0;
box_length = 40;
box_height = 15;
//
car_cl = 0.05;
car_offset_length = 3;
car_offset_height = 0.2;
car_length = 5.5;
car_height = 1.4;
car_diffuser_length = 1.2;
car_diffuser_height = 0.3;
//
// geometry [box, car]
//
Point(1) = {0, 0, 0, box_cl};
Point(2) = {car_offset_length-car_length*0.2, 0, 0, car_cl};
Point(3) = {car_offset_length+car_length*1.2, 0, 0, car_cl};
Point(4) = {box_length, 0, 0, box_cl};
Point(5) = {box_length, box_height, 0, box_cl};
Point(6) = {0, box_height, 0, box_cl};
//
Point(7) = {car_offset_length, car_offset_height, 0, car_cl};
Point(8) = {car_offset_length+car_length-car_diffuser_length, car_offset_height, 0, car_cl};
Point(9) = {car_offset_length+car_length, car_offset_height+car_diffuser_height, 0, car_cl};
Point(10) = {car_offset_length+car_length, car_offset_height+0.7, 0, car_cl};
Point(11) = {car_offset_length+car_length*0.8, car_offset_height+car_height*0.7, 0, car_cl};
Point(12) = {car_offset_length+car_length*0.7, car_offset_height+car_height, 0, car_cl}; // roof
Point(13) = {car_offset_length+car_length*0.5, car_offset_height+car_height, 0, car_cl}; // roof
Point(14) = {car_offset_length+car_length*0.3, car_offset_height+car_height*0.5, 0, car_cl};
Point(15) = {car_offset_length + 0.1, car_offset_height+car_height*0.3, 0, car_cl};
//
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};
//
Line(7) = {7, 8};
Line(8) = {8, 9};
// Line(9) = {9, 10};
// Line(10) = {10, 11};
// Line(11) = {11, 12};
// Line(12) = {12, 13};
// Line(13) = {13, 14};
// Line(14) = {14, 15};
// Line(15) = {15, 7};
//
Spline(16) = {9, 10, 11};
Spline(17) = {11, 12, 13, 14};
Spline(18) = {14, 15, 7};
//
Curve Loop(1) = {1, 2, 3, 4, 5, 6};
//
// Curve Loop(2) = {12, 13, 14, 15, 7, 8, 9, 10, 11};
Curve Loop(2) = {17, 18, 7, 8, 16};
//
Plane Surface(1) = {1, 2};
//
// boundary conditions [box, car]
//
Physical Line(11) = {6};
Physical Line(21) = {4};
Physical Line(31) = {1,2,3};
Physical Line(41) = {5};
//
// Physical Line(32) = {7,8,9,10,11,12,13,14,15};
Physical Line(32) = {7,8,16,17,18};
//
// domain conditions [fluid]
//
Physical Surface(100) = {1};
//

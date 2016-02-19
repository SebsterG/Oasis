// Gmsh project created on Fri Feb 19 10:15:08 2016
Point(1) = {-1, 0, 0, 0.0001};
Point(2) = {1, 0, 0, 0.0001};
Point(3) = {0.0, -4, 0, 0.04};
Point(4) = {0.5, -2.0, 0, 0.04};
Point(5) = {-0.5, -2.0, 0, 0.04};
Line(1) = {1, 5};
Line(2) = {5, 3};
Line(3) = {3, 4};
Line(4) = {4, 2};
Line(5) = {2, 1};
Point(6) = {0.0, 0.0, 0, 0.01};
Line(6) = {1, 6};
Line(7) = {6, 2};
Line Loop(8) = {1, 2, 3, 4, -7, -6};
Plane Surface(9) = {8};

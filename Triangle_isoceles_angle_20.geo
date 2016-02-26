// Gmsh project created on Fri Feb 19 13:45:56 2016
Point(1) = {0.0, 1.0, 0, 1
.0};
Point(2) = {1.0, 1.0, 0, 1.0};
Point(3) = {1.0, 0.0, 0, 1.0};
Point(4) = {1, 0.5, -0, 1.0};
Point(5) = {0.5, 1.0, -0, 1.0};
Point(6) = {0.5, -0.5, -0, 1.0};
Delete {
  Point{6};
}
Point(6) = {0.5, 0.5, -0, 1.0};
Line(1) = {1, 6};
Line(2) = {3, 3};
Line(3) = {6, 3};
Line(4) = {3, 4};
Line(5) = {4, 2};
Line(6) = {2, 5};
Line(7) = {5, 1};
Line Loop(8) = {7, 1, 3, 4, 5, 6};
Plane Surface(9) = {8};

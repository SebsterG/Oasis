cl__1 = 0.00098;
cl__2 = 0.0098;
Point(1) = {0, -2, 0, cl__2};
Point(2) = {1.732050807568877, 1, 0, cl__1};
Point(3) = {-1.732050807568877, 1, 0, cl__1};
Point(4) = {0, 1, 0, cl__2};
Point(6) = {0.8660254037844386, -0.5, 0, cl__2};
Point(8) = {-0.8660254037844386, -0.5, 0, cl__2};
Line(3) = {3, 8};
Line(4) = {8, 1};
Line(5) = {1, 6};
Line(6) = {6, 2};
Line(7) = {2, 4};
Line(8) = {4, 3};
Line Loop(10) = {3, 4, 5, 6, 7, 8};
Plane Surface(10) = {10};

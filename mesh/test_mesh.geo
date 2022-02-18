dx1 = 0.05;
dx2 = 0.01;
dx3 = 0.025;

circ_rad = 0.25;

Point(1) = {0, 0, 0, dx1};
Point(2) = {2.25-circ_rad, 0, 0, dx2};
Point(3) = {2.25+circ_rad, 0, 0, dx2};
Point(4) = {3.5, 0, 0, dx1};
Point(5) = {3.5, 2, 0, dx1};
Point(6) = {2.25+circ_rad, 2, 0, dx2};
Point(7) = {2.25-circ_rad, 2, 0, dx2};
Point(8) = {0, 2, 0, dx1};

Point(9) = {1.25, 0.5, 0, dx2};
Point(10) = {1.25, 0.5-circ_rad, 0, dx2};
Point(11) = {1.25-circ_rad, 0.5, 0 ,dx2};
Point(12) = {1.25, 0.5+circ_rad, 0, dx2};
Point(13) = {1.25+circ_rad, 0.5, 0, dx2};

Point(14) = {1.25, 1.5, 0, dx2};
Point(15) = {1.25, 1.5-circ_rad, 0, dx2};
Point(16) = {1.25-circ_rad, 1.5, 0 ,dx2};
Point(17) = {1.25, 1.5+circ_rad, 0, dx2};
Point(18) = {1.25+circ_rad, 1.5, 0, dx2};

Point(19) = {2.25, 1, 0, dx2};
Point(20) = {2.25, 1-circ_rad, 0, dx2};
Point(21) = {2.25-circ_rad, 1, 0 ,dx2};
Point(22) = {2.25, 1+circ_rad, 0, dx2};
Point(23) = {2.25+circ_rad, 1, 0, dx2};

Point(24) = {2.25, 0, 0, dx2};
Point(25) = {2.25-circ_rad, 0, 0 ,dx2};
Point(26) = {2.25, circ_rad, 0, dx2};
Point(27) = {2.25+circ_rad, 0, 0, dx2};

Point(28) = {2.25, 2, 0, dx2};
Point(29) = {2.25-circ_rad, 2, 0 ,dx2};
Point(30) = {2.25, 2-circ_rad, 0, dx2};
Point(31) = {2.25+circ_rad, 2, 0, dx2};
//+
Point(32) = {1.75, 0, 0, dx2};
//+
Point(33) = {1.75, 2, 0, dx2};
//+
Point(34) = {2.75, 0, 0, dx1};
//+
Point(35) = {2.75, 2, 0, dx3};
//+
Line(2) = {32, 2};
//+
Circle(3) = {2, 24, 26};
//+
Circle(4) = {26, 24, 3};
//+
Line(5) = {3, 34};
//+
Line(6) = {34, 4};
//+
Line(7) = {4, 5};
//+
Line(8) = {5, 35};
//+
Line(9) = {35, 6};
//+
Circle(10) = {6, 28, 30};
//+
Circle(11) = {30, 28, 7};
//+

//+
Circle(17) = {20, 19, 23};
//+
Circle(18) = {23, 19, 22};
//+
Circle(19) = {22, 19, 21};
//+
Circle(20) = {21, 19, 20};
//+
Circle(21) = {13, 9, 12};
//+
Circle(22) = {12, 9, 11};
//+
Circle(23) = {11, 9, 10};
//+
Circle(24) = {10, 9, 13};
//+
Circle(25) = {15, 14, 18};
//+
Circle(26) = {18, 14, 17};
//+
Circle(27) = {17, 14, 16};
//+
Circle(28) = {16, 14, 15};

//+
//+
//+
//+
Point(38) = {1.25, 0, 0, dx1};
//+
Point(39) = {1.25, 2, 0, dx3};
//+
Line(29) = {1, 38};
//+
Line(30) = {38, 32};
//+
Line(31) = {7, 33};
//+
Line(32) = {33, 39};
//+
Line(33) = {39, 8};
//+
Line(34) = {8, 1};
//+
Line(35) = {38, 10};
//+
Line(36) = {12, 15};
//+
Line(37) = {17, 39};
//+
Line(38) = {30, 22};
//+
Line(39) = {20, 26};
//+
Line Loop(1) = {29, 35, -23, -22, 36, -28, -27, 37, 33, 34};
//+
Plane Surface(1) = {1};
//+
Line Loop(2) = {32, -37, -26, -25, -36, -21, -24, -35, 30, 2, 3, -39, -20, -19, -38, 11, 31};
//+
Plane Surface(2) = {2};
//+
Line Loop(3) = {7, 8, 9, 10, 38, -18, -17, 39, 4, 5, 6};
//+
Plane Surface(3) = {3};
//+
//+

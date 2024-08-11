D = 0.2;
L = 22.0*D;
H = 4.1*D;

hCoarse = 0.05;
hFine   = 0.05*hCoarse;
hFin = 0.01;



// Outer rectangle points
Point(1)  = {0.0, 0.0, 0.0, hCoarse};
Point(2)  = {0.0, H,   0.0, hCoarse};
Point(3)  = {L,   H,   0.0, hCoarse};
Point(4)  = {L,   0.0, 0.0, hCoarse};



// Left circle points
Point(5)  = {2.0*D, 2.0*D, 0.0, hFine};
Point(6)  = {1.5*D, 2.0*D, 0.0, hFine};
Point(7)  = {2.5*D, 2.0*D, 0.0, hFine};



// Outer rectangle lines (both inlets and top lines)
Line(2)   = {3, 4};
Line(4)   = {1, 2};



// Left circle defined
Circle(5) = {6, 5, 7};
Circle(6) = {7, 5, 6};



// Funnel points bottom defined
Point(8) = {2.0, 0, 0, hCoarse};
//+
Point(9) = {2.4, 0, 0, hCoarse};
//+
Point(10) = {2.0, -0.82, 0, hCoarse};
//+
Point(11) = {2.4, -0.82, 0, hCoarse};



// Funnel points top defined
Point(17) = {2.0, 0.82, 0, hCoarse};
//+
Point(18) = {2.4, 0.82, 0, hCoarse};
//+
Point(19) = {2.0, 1.64, 0, hCoarse};
//+
Point(20) = {2.4, 1.64, 0, hCoarse};



// Funnel lines and lines to connect to outer rectangle
Line(7) = {1, 8};
//+
Line(8) = {8, 10};
//+
Line(9) = {10, 11};
//+
Line(10) = {11, 9};
//+
Line(11) = {9, 4};



// Right Cylinder Points
Point(12) = {4.0, 0.4, 0, hFine};
//+
Point(13) = {3.9, 0.4, 0, hFine};
//+
Point(14) = {4.1, 0.4, 0, hFine};



// Right Cylinder Mesh Cutout
Circle(15) = {13, 12, 14};
Circle(16) = {14, 12, 13};



// Top Funnel Line Connections
Line(17) = {2, 17};
//+
Line(18) = {17, 19};
//+
Line(19) = {19, 20};
//+
Line(20) = {20, 18};
//+
Line(21) = {18, 3};

// Fin Points
Point(21) = {0.51, 0.45, 0, hFin};
//+
Point(22) = {0.51, 0.35, 0, hFin};
//+
Point(23) = {1.5, 0.35, 0, hFin};
//+
Point(24) = {1.5, 0.45, 0, hFin};
//+
Point(25) = {3.89, 0.45, 0, hFin};
//+
Point(26) = {3.89, 0.35, 0, hFin};
//+
Point(27) = {2.9, 0.35, 0, hFin};
//+
Point(28) = {2.9, 0.45, 0, hFin};



Line(22) = {22, 23};
//+
Line(23) = {23, 24};
//+
Line(24) = {24, 21};
//+
Line(25) = {21, 22};
//+
Line(26) = {28, 25};
//+
Line(27) = {25, 26};
//+
Line(28) = {26, 27};
//+
Line(29) = {27, 28};
//+
Curve Loop(1) = {17, 18, 19, 20, 21, 2, -11, -10, -9, -8, -7, 4};
//+
Curve Loop(2) = {6, 5};
//+
Curve Loop(3) = {24, 25, 22, 23};
//+
Curve Loop(4) = {26, 27, 28, 29};
//+
Curve Loop(5) = {16, 15};
//+
Plane Surface(1) = {1, 2, 3, 4, 5};
//+

// Physical Curve IDs
Physical Curve(30) = {4}; // Inlet Left
//+
Physical Curve(31) = {17}; // Top Left
//+
Physical Curve(32) = {18}; // Funnel Top Left Wall
//+
Physical Curve(33) = {19}; // Outlet Top
//+
Physical Curve(34) = {20}; // Funnel Top Right Wall
//+
Physical Curve(35) = {21}; // Top Right
//+ 
Physical Curve(36) = {2}; // Inlet Right
//+
Physical Curve(37) = {11}; // Bottom Right
//+
Physical Curve(38) = {10}; // Funnel Bottom Right Wall
//+
Physical Curve(39) = {9}; // Outlet Bottom
//+
Physical Curve(40) = {8}; // Funnel Bottom Left Wall
//+
Physical Curve(41) = {7}; // Bottom Left
//+
Physical Curve(42) = {5, 6}; // Cylinder Left
//+
Physical Curve(43) = {22, 23, 24, 25}; // Fin Left
//+
Physical Curve(44) = {26, 29, 28, 27}; // Fin Right
//+
Physical Curve(46) = {16, 15}; // Cylinder Right
//+
Physical Surface(47) = {1};

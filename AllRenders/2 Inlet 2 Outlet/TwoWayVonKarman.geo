D = 0.2;
L = 22.0*D;
H = 4.1*D;

hCoarse = 0.05;
hFine   = 0.05*hCoarse;



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
//+
Curve Loop(1) = {4, 17, 18, 19, 20, 21, 2, -11, -10, -9, -8, -7};
//+
Curve Loop(2) = {5, 6};
//+
Curve Loop(3) = {15, 16};
//+
Plane Surface(1) = {1, 2, 3};
//+
Physical Curve(22) = {4}; // Inlet Left
//+
Physical Curve(23) = {17}; // Top Left
//+
Physical Curve(24) = {18}; //Funnel Top Left Wall
//+
Physical Curve(25) = {19}; // Outlet Top
//+
Physical Curve(26) = {20}; // Funnel Top Right Wall
//+
Physical Curve(27) = {21}; // Top Right
//+
Physical Curve(28) = {2}; // Inlet Right
//+
Physical Curve(29) = {11}; // Bottom Right
//+
Physical Curve(30) = {10}; // Funnel Bottom Right Wall
//+
Physical Curve(31) = {9}; // Outlet Bottom
//+
Physical Curve(32) = {8}; // Funnel Bottom Left Wall
//+
Physical Curve(33) = {7}; // Bottom Left
//+
Physical Curve(34) = {5, 6}; // Circle Left
//+
Physical Curve(35) = {15, 16}; // Circle Right
//+
Physical Surface(36) = {1};

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
Line(1)   = {2, 3};
Line(2)   = {3, 4};
Line(4)   = {1, 2};



// Left circle defined
Circle(5) = {6, 5, 7};
Circle(6) = {7, 5, 6};



// Funnel points defined
Point(8) = {1.8, 0, 0, hFine};
//+
Point(9) = {2.6, 0, 0, hFine};
//+
Point(10) = {1.8, -0.82, 0, hCoarse};
//+
Point(11) = {2.6, -0.82, 0, hCoarse};



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



// Defining outer loop and surface
Curve Loop(1) = {1, 2, -11, -10, -9, -8, -7, 4};
//+
Curve Loop(2) = {6, 5};
//+
Curve Loop(3) = {16, 15};
//+
Plane Surface(1) = {1, 2, 3};



// Physical Curves (aka IDs for code)
Physical Curve(17) = {4}; // Inlet left
//+
Physical Curve(18) = {1}; // Top
//+
Physical Curve(19) = {2}; // Inlet Right
//+
Physical Curve(20) = {11}; // Bottom Right
//+
Physical Curve(21) = {10}; // Funnel Right
//+
Physical Curve(22) = {9}; // Outlet (Bottom)
//+
Physical Curve(23) = {8}; //Funnel Left
//+
Physical Curve(24) = {7}; // Bottom Left
//+
Physical Curve(26) = {5, 6}; // Circle Left
//+
Physical Curve(27) = {16, 15}; // Circle Right
//+
Physical Surface(25) = {1};

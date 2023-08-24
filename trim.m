% prints the trim conditions

%% params

m = 1.5;g = 9.81; S = 0.288; b = 1.2; cb = 0.24; 
c_d0 = 0.0065; k = 0.06366; cL_0 = 0.3465; cL_a = 4.198; cL_q = 2.2423; cL_de = 0.3862;
cm_0 = 0.02; cm_a = -0.5456; cm_q = -5.1639; cm_de = -0.8895;
cy_0 = 0; cy_b = -0.1623; cy_p = 0; cy_r = 0.1495; cy_dr = 0.02575;
cl_0 = 0; cl_b = -0.00217; cl_p = 0; cl_r = 0; cl_dr = 0.0027; cl_da = 0;
cn_0 = 0; cn_b = 0.148255; cn_p = 0; cn_da = 0; cn_r = -0.0688; cn_dr = -0.0593;
%%
z = 0;% m
rho0 = 1.225; % kg/m^3
H0 = 10400; % m (scale height)

rho = rho0*exp(-z/H0);

Vinf = 21.23; % m/s---------------------

W = m*g;
Q = 1/2*rho*Vinf^2;

cLtrim = W / (Q*S)
cDtrim = c_d0 + k*cLtrim^2;
thrust = Q*S*cDtrim

res = ([cL_a cL_de; cm_a cm_de])\[cLtrim-cL_0; -cm_0];

alpha_trim = res(1)
dele_trim = res(2)

u_trim = Vinf*cos(alpha_trim)
w_trim = Vinf*sin(alpha_trim)
winstyle = 'docked';
% winstyle = 'normal';

set(0,'DefaultFigureWindowStyle',winstyle)
set(0,'defaultaxesfontsize',18)
set(0,'defaultaxesfontname','Times New Roman')
% set(0,'defaultfigurecolor',[1 1 1])

% clear VARIABLES;
clear
global spatialFactor;
global c_eps_0 c_mu_0 c_c c_eta_0
global simulationStopTimes;
global AsymForcing
global dels
global SurfHxLeft SurfHyLeft SurfEzLeft SurfHxRight SurfHyRight SurfEzRight

%Jinseng Vanderkloot 101031534

% Question 2.a: This simulates a 2D FDFT for a traveling electro magnetic 
% wave from a source which trvels into a different inclusion (different 
% permiativity) using Yee cells calculation. 

% Question 2.b.i: 
%Commenting second line with ep{1} removes inclusion and the wave just
%travels left to right (waves expanding outwards due to open area) 
%with no barriers. 

% Question 2.b.ii: 
% BC defines the boundary conditions for the area as well as defines the
% source information 

% Question 2.b.iii:
% bcf1g.s(1) is setting up the wave source parameters like location, type
% and magnitude. In PlaneWaveBC.m there are two types "p" and "s" which
% change the inputted wave type. Attempted to change Phi, omega and betap
% which changes wave source parameters for s type but dont know what values
% are appropriate in input, didnt see much change. The value of s=1 or 100  
% changed refelction or scattering of wave I think since it bounces back 
% between source and inclusion in EZ and Hy. 

%Question 2.b.iv:
% When changing bc{1}.xp to 'e', the wave is found bouncing back
% all around area = no disipation on the boundary (reflectiong) and wave 
% continues to travel and scatter within region with no escape. 

%Question 3 a: added more inclusions 
%Question 3 b: Change st to -0.05 = Grading is more obvious 
%Quesiton 3 c: Double frequnecy = more reflection, half frequency = more
%wave bouncing within inclusion 

%Question 4: Added 2 more inclusiong to 'trap' wave making Q3C more obvious
%about the refelction and wave traping by inclusions when increasing and 
%decreasing freqeuncy.   

%Inital Parameters 
dels = 0.75;
spatialFactor = 1;

c_c = 299792458;                  % speed of light
c_eps_0 = 8.8542149e-12;          % vacuum permittivity
c_mu_0 = 1.2566370614e-6;         % vacuum permeability
c_eta_0 = sqrt(c_mu_0/c_eps_0);


tSim = 200e-15
f = 2* 230e12;
lambda = c_c/f;

xMax{1} = 20e-6;
nx{1} = 200;
ny{1} = 0.75*nx{1};

%All region setup 
Reg.n = 1;

mu{1} = ones(nx{1},ny{1})*c_mu_0;

%Inclusion Setup 
epi{1} = ones(nx{1},ny{1})*c_eps_0;
epi{1}(150:180,55:95)= c_eps_0*11.3;
epi{1}(60:125,95:110)= c_eps_0*11.3;
epi{1}(60:125,45:55)= c_eps_0*11.3;

sigma{1} = zeros(nx{1},ny{1});
sigmaH{1} = zeros(nx{1},ny{1});

dx = xMax{1}/nx{1};
dt = 0.25*dx/c_c;
nSteps = round(tSim/dt*2);
yMax = ny{1}*dx;
nsteps_lamda = lambda/dx

%Plotting Setup 
movie = 1;
Plot.off = 0;
Plot.pl = 0;
Plot.ori = '13';
Plot.N = 100;
Plot.MaxEz = 1.1;
Plot.MaxH = Plot.MaxEz/c_eta_0;
Plot.pv = [0 0 90];
Plot.reglim = [0 xMax{1} 0 yMax];


%Boundary Condition/Wave Source setup 
bc{1}.NumS = 1;
bc{1}.s(1).xpos = nx{1}/(4) + 1;
bc{1}.s(1).type = 'ss';
bc{1}.s(1).fct = @PlaneWaveBC;
%mag = -1/c_eta_0;
mag = 1;
phi = 1;
omega = f*2*pi;
betap = 1;
t0 = 30e-15;
st = -0.05;
s = 0;
y0 = yMax/2;
sty = 1.5*lambda;
bc{1}.s(1).paras = {mag,phi,omega,betap,t0,st,s,y0,sty,'s'};

Plot.y0 = round(y0/dx);

bc{1}.xm.type = 'a';
bc{1}.xp.type = 'a';
bc{1}.ym.type = 'a';
bc{1}.yp.type = 'a';

pml.width = 20 * spatialFactor;
pml.m = 3.5;

Reg.n  = 1;
Reg.xoff{1} = 0;
Reg.yoff{1} = 0;

RunYeeReg







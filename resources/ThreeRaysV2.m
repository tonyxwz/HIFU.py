% computation of power (in watt) through triagle, perpendicular on triple
% ray, version with simplified computation, distance to triangle is 
% not needed anymore

clear all
close all

r1=[1,2,3]; % arbitrary position of the transducer
n1=[1,0,1]; % arbitrary normal on transducer = direction of beam
n1=n1/norm(n1); % normalized n1


v1=[1.4,0.4,0.3]; % carrier ray direction
v1=v1/norm(v1); % normalized v1
 
% check if v1 does not fire backwards:
cosphi=v1*n1';
phi=acos(cosphi);
phi_degrees=phi/pi*180;
fprintf('phi= %6.2f deg\n',phi_degrees);
if cosphi<0
    error('ray direction is backwards');
end;   

% generate two other rays
b2=cross(v1,[1,0,0]); % [1,0,0] =arbitrary vector, not parallel to v1
b3=cross(v1,b2);
b2=b2/norm(b2);
b3=b3/norm(b3);
% b2 and b2 are perpendicular to each other and to v1
gamma=0.05;  % angle between v1 and v2 and between v1 and v3
fprintf('angle between rays and carrier ray: gamma = %8.3e deg\n', gamma/pi*180)
% gamma must be much smaller when using in ray tracer,
% but then you cannot see anything anymore on the plots
v2=cos(gamma)*v1 + sin(gamma)*b2;
v3=cos(gamma)*v1 + sin(gamma)*b3;
% v1, v2  and v3 form the triple ray

% plot triple ray
lambda=5; % length to plot
u1=r1+lambda*v1;
u2=r1+lambda*v2;
u3=r1+lambda*v3;
figure
plot3(r1(1),r1(2),r1(3),'*r');% plot position of tranducer 
hold on;

% plot normal on transducer over distance beta
beta=1;
plot3([r1(1),r1(1)+beta*n1(1)],[r1(2),r1(2)+beta*n1(2)],[r1(3),r1(3)+beta*n1(3)],'-m')

% plot the three rays, carrier ray in red:
xyz=[r1;u1];
xx1=xyz(:,1); yy1=xyz(:,2); zz1=xyz(:,3);
xyz=[r1;u2];
xx2=xyz(:,1); yy2=xyz(:,2); zz2=xyz(:,3);
xyz=[r1;u3];
xx3=xyz(:,1); yy3=xyz(:,2); zz3=xyz(:,3);
plot3(xx1,yy1,zz1,'-r',xx2,yy2,zz2,'-b',xx3,yy3,zz3,'-g');

% take a point on distance 1 on the carrier ray
% ditance P1 to transducer
% P1=r1+v1; % point on carrier ray with distance lambda to transducer element
% compute points P2 and P3 where the other two rays pass through the plane that is
% perpendicular to the carrier ray in P1
% find P2=r1+ mu2*v2 such that P2-P1 is perpendicular to v1
% so (P2-P1)*v1'=0, hence (r1+mu2*v1-P1)*v1'=0 
% mu2= ( (P1-r1)*v1' ) /( v2*v1');
% P2= r1+mu2*v2;
% same for third ray:
% mu3= ( (P1-r1)*v1' ) /( v3*v1');
% P3= r1+mu3*v3;

% More clever computation:
% the angle between v1 and v2 is gamma, hence the length of P1-P2 is
d12= tan(gamma);
% also the angle between v1 and v3 is gamma, so the length of P1-P3 is
d13= tan(gamma);
% the cos of the angle between v2 and v3 is v2*v3', hence
d23= tan(acos(v2*v3'));

s=(d12+d13+d23)/2;
% Heron's formula
Atriangle=sqrt(s*(s-d12)*(s-d13)*(s-d23));
fprintf('the area of the triangle (at distance 1) is %10.5e \n',Atriangle);

% compute intensity using far field approx,\
% parameters may be changed
RequiredTotalPower=1; % Watt
fprintf('RequiredTotalPower: %5.1f watt\n',RequiredTotalPower);
Radius = 0.002;% m
Atrd=pi*Radius^2; % transducer area
velocity=15.37*100; % speed of sound in gel, in: m/s
freq = 1.2e6; % Hz;
k= 2*pi*freq/velocity; % wave number
density=1040; % gel, density in kg/m^3
Zacous = density*velocity;
fprintf('Trans elementRadius= %7.4f m, wave length lambda= %9.6f m\n',Radius,2*pi/k);
S0=ComputeS0(Radius,k,Zacous,RequiredTotalPower);

% far field approx in P1 (distance 1 to transducer element
% pressure p(1, phi) = Atrd S0 exp(i k 1)/(2 pi 1)  D(theta) with:
% Atrd= transducer area=pi Radius^2 
% D(theta)= 2 J1(Radius*k sin(phi)) / (Radius*k*sin(phi))
% J1=Besselfunction, n=1
% intensity in P1: I = |p|^2/(2*Z)
ka=k*Radius;
D=(2*besselj(1,ka*sin(phi))/(ka*sin(phi)));
I= Atrd^2*S0^2/(8*pi^2*Zacous)*D^2; % intensity at P1, distance 1

% POWER THROUGH TRIANGLE:
Pow=I*Atriangle; % power going through the triangle
fprintf('power through the triangle %10.5e\n',Pow);

xlabel('x'); ylabel('y'); zlabel('z');
disp('done');

ComputeS0(3.3e-3, k,Zacous, 1)
getS0(Radius, freq, density, velocity, RequiredTotalPower)

function S0=ComputeS0(Radius,k,Z,RequiredPower)
   % R=transd element radius, k=wave number=2 pi/lambda
   % Z= acoustic impedance, 
   % computes S0 such that flat transducer element generates RequiredPower
   % pressure field of this element (far field approx):
   % p(r, theta) = A S0 exp(i k r)/(2 pi r)  D(theta) with:
   % A= area=pi Radius^2
   % D(theta)= 2 J1(Radius*k sin(theta)) / (Radius*k*sin(theta))
   % D = directivity = determines pressure field in direction theta
   % J1=Besselfunction, n=1
   ka=Radius*k;
   A= pi*Radius^2;
   f=@(theta) sin(theta).*(2*besselj(1,ka*sin(theta))./(ka*sin(theta))).^2;
   Int=quad(f,0.00000001,pi/2); 
   S0= sqrt(4*pi*Z*RequiredPower/Int)/A;
   % done
end

function s0 = getS0(R, freq, density, speed, RequiredPower)
    ka=Radius*
    fx = @(theta) sin(theta).*abs(besselj(1,ka*sin(theta))./sin(theta)).^2;
    maxflux = integral(fx, 0, pi/2);
    s0 = RequiredPower*4*(freq^2)*pi*density/...
		 (speed * R^2 * maxflux);
end




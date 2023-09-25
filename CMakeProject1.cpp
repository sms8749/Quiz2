clc;
clear;
%% Main Code - Problem 1
%Need to check to see if i get same values as from project report 
%also r == r_2
% correct eci to ecef frame for 
GMST_1 = 225.679;
r = [1138.118*cosd(GMST_1) + 5354.35*sind(GMST_1); 1138.11*sind(GMST_1)-5354.36*cosd(GMST_1); 4758.45];
v = [7.775*-0.5204; 7.775*-0.8537; 7.775*0.0207];
mu = (3.986 * 10^5);

a_2 = 8059; %km
e_2 = 0.1; %No units
I_2 = 41; %deg
RAAN_2 = 60; %deg
AOP_2 = 90.0; %deg
f_2 = 0; %deg

[a, e, I, RAAN, AOP, f] = RV2OE(r,v,mu);
fprintf ('a = %0.3f km \n', a);
fprintf ('e = %0.3f \n', e);
fprintf ('I = %0.3f deg\n', I);
fprintf ('RAAN = %0.3f deg \n', RAAN);
fprintf ('AOP = %0.3f deg \n', AOP);
fprintf ('f = %0.3f deg\n', f);

[r_2,v_2] = OE2RV (a_2, e_2, I_2, RAAN_2, AOP_2, f_2, mu);
fprintf ('r in ECI frame = [%0.3f, %0.3f, %0.3f] km\n', r_2);
fprintf ('v in ECI frame = [%0.3f, %0.3f, %0.3f] km/s\n', v_2);

r_eci = r_2';
v_eci = v_2';
omega = 41; %deg
t = 0; %time
GMST = 225.679; %deg MUST FIX

[r_ecef, v_ecef] = ECI2ECEF(r_eci, v_eci, omega, t, GMST, mu);
fprintf ('r in ECEF frame = [%0.3f, %0.3f, %0.3f] km\n', r_eci);
fprintf ('v in ECEF frame = [%0.3f, %0.3f, %0.3f] km/s\n', v_eci);

%% Main Code - Problem 2
% 1. Propage orbit of satellite in the ECI frame using ODE45
% 1. Plot orbit of the satellite in the ECI frame

% 2. Plot same orbit in ECEF frame in:
% 3D, X-Y plane view, Y-Z plane view, X-Z plane view --> comment on results

% 3. Plot ground track for 24 orbits and location of SC.
% How many times does it pass over SC? why?
% Does this change with change in eccentricity, inclination, and RAAN? why?
% Demonstrate by using another orbit with sufficiently diff inclination and
% eccentircty. can you explain the ground track?

% 4. Repeat (1.) using F&G functions --> this is analytic propagation

% 5. Compare errors in r and v between numerical and analytic propagation.
% Comment on mag of errors. expected? which is more accurate?

% 6. Compute orbital elements from r and v numerically propagated
% expected? what do you observe?


%% Functions - Part 1
% 

% Function 1: RV2OE
% Input: R + V in ECI frame, mu 
% Output: Orbital Elements
function [a, e, I, RAAN, AOP, f] = RV2OE(r, v, mu)

mag_r = sqrt(sum(r.^2)); %mag of r
E = (dot(v, v)/2) - (mu/mag_r); % E formula
a = -mu/(2*E); %Find alpha

h_cross = cross(r, v); %h vector
ev = (cross(v, h_cross)/mu) - (r/mag_r); %e vector
e = sqrt(sum(ev.^2)); %e mag
    
i = [1;0;0];
j = [0;1;0];
k = [0;0;1];

mag_h = sqrt(sum(h_cross.^2)); %mag of h
I = acosd(dot(k,h_cross)/(mag_h)); %i 

h = cross(r, v); %h vector
k_h = sqrt(sum(cross(k,h).^2)); %K cross h vector
n = cross(k,h)/ k_h; %n formula
check = dot(i,n); %check for n dot I is less than 0
if check < 0 % less than zero add pi
    RAAN = atand(dot(j,n)/(dot(i,n))) + pi;
else % greater than zero, nothing else
    RAAN = atand(sin(dot(j,n))/(dot(i,n)));
end    

h_cross = cross(r, v); %h vecoter
e_vec = (cross(v, h_cross)/mu) - (r/mag_r); % e vector
e = sqrt(sum(e_vec.^2)); % mag of e
h = cross(r, v); % h vector
Vk_h = cross(k,h); % K cross h (vector)
k_h = sqrt(sum(Vk_h.^2)); % mag of K cross h
n = cross(k,h)/ k_h; % n formula
check = dot(e_vec, k); %if n dot K is less than 0
if check < 0 % less than zero negative sign
    AOP = -acosd(dot(e_vec, n)/ e);
else % greater than zero, leave alone
    AOP = acosd(dot(e_vec, n)/ e);
end


    if e == 0
        f =0;
    else
        f = acos(dot(ev, r) / e * norm(r));
        if dot(r,v) < 0
            f = 2 * pi - f;
        end
    end

end


% Function 2: OE2RV (must have)
% Input: Orbital Elements and mu
% Output: R + V in ECI frame
function [r_2,v_2] = OE2RV (a_2, e_2, I_2, RAAN_2, AOP_2, f_2, mu)

h = sqrt(mu * a_2 * (1-e_2^2));

r_ecef = (h^2/mu) * (1/(1+ e_2 * cosd(f_2))) * [cosd(f_2); sind(f_2); 0]; %store r ecef as a vector

dcm_ecef_to_eci = [ cosd(RAAN_2)*sind(AOP_2)-sind(RAAN_2)*sind(AOP_2)*cosd(I_2),   sind(RAAN_2)*cosd(AOP_2)+cosd(RAAN_2)*sind(AOP_2)*cosd(I_2),   sind(RAAN_2)*sind(I_2);
                   -cosd(RAAN_2)*sind(AOP_2)-sind(AOP_2)*cosd(AOP_2)*cosd(I_2),   -sind(RAAN_2)*sind(AOP_2) + cosd(RAAN_2)*cosd(AOP_2)*cosd(I_2), cosd(AOP_2)*sind(I_2);
                                   sind(AOP_2)*sind(I_2),                            -cosd(RAAN_2)*cosd(I_2),                                  cosd(I_2)];
r_2 = dcm_ecef_to_eci * r_ecef;

v_mag_2 = sqrt(mu*(2/norm(r_2) - 1/a_2));
v_ecef = (mu/h) * [-sind(f_2); e_2+cosd(f_2);0];
v_2 = dcm_ecef_to_eci * v_ecef;

r_2 = r_2';
v_2 = v_2';

    end


% Function 3: ECI2ECEF (must have)
% Input: R + V in ECI Frame, GMST, mu, time, omega ( w )
function[r_ecef, v_ecef] = ECI2ECEF(r_eci, v_eci, omega, t, GMST, mu)

dcm_eci_ecef_GMST = [cosd(GMST+(omega*t)),  sind(GMST+(omega*t)), 0;
                     -sind(GMST+(omega*t)), cosd(GMST+(omega*t)), 0;
                     0,                     0,                    1];
r_ecef = dcm_eci_ecef_GMST * r_eci;
r_ecef = r_ecef';
%v_ecef_rotation = cross([0;1;0], r_eci);
v_ecef = (dcm_eci_ecef_GMST .* [0;1;0]) * v_eci;
v_ecef= v_ecef';

    end

%% Functions for Problem 2
% Function 4: FG to find F & G Functions of Lagrange
% Input: inital position, velocity vector, and time elapsed
% Output: final position and velocity vector
function [r, v] = FGFunc(r0,v0,dt,mu)


    end


% Function 5: codes equations of motion for 2 body problem to use with ode45
% Input:
% Output:
function [dx] = TwoBP(t, x, mu)


    end


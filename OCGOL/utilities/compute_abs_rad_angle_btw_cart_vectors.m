function [theta_abs] = compute_abs_rad_angle_btw_cart_vectors(u,v)

%cartesian inputs
%2 x,y vectors

%% Define vectors to calculate angle between

uCar = u;
vCar = v;


%% Calulate angle between 2 vectors (radians)
theta = atan2(uCar(1)*vCar(2) - uCar(2)*vCar(1), uCar(1)*vCar(1) + uCar(2)*vCar(2));
%get absolute distance (smallest angle)
theta_abs = abs(theta);


end


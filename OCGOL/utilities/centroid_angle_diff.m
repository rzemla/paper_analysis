function [theta] = centroid_angle_diff(car_vectors, rewardLoc)
%returns the angle difference between reward location and ts centroid of
%each ROI

    %convert reward location to complex vector and then cartesian
    %get a track length normalized output from Behavior struct
    u_rew = exp(1i*deg2rad((rewardLoc/196)*360));
    uCar_rew = [real(u_rew), imag(u_rew)];
    
    %plot reward location
    figure;
    compass(u_rew)
    
    %for each ROI get difference between ts vector and reward vector
    for ROI = 1:size(car_vectors,2)
        uCar{1} = car_vectors{ROI};
        uCar{2} = uCar_rew;
        
        %get smallest angle between vectors
        theta(ROI) = atan2(uCar{1}(1)*uCar{2}(2) - uCar{1}(2)*uCar{2}(1), uCar{1}(1)*uCar{2}(1) + uCar{1}(2)*uCar{2}(2));
        theta(ROI) = abs(theta(ROI));
    end

    %% Verify a random ROI with a compass plot - TODO
    
end


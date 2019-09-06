function [theta_abs] = centroid_angle_diff_btw_ROIs_partial(car_vectors, car_vectors_2,pf_vector_max)
%returns the angle difference between reward location and ts centroid of
%each ROI

    %convert reward location to complex vector and then cartesian
    %get a track length normalized output from Behavior struct
    %u_rew = exp(1i*deg2rad((rewardLoc/196)*360));
    %uCar_rew = [real(u_rew), imag(u_rew)];
    
    %plot reward location
    figure;
    %compass(u_rew)
    
    
    %for each ROI get difference between ts vector and reward vector
    for ROI = 1:size(car_vectors,2)
        %         uCar{1} = car_vectors{ROI};
        %         uCar{2} = uCar_rew;
        %uCar{1} = car_vectors{ROI};
        %uCar{2} = car_vectors_2{ROI};
        
        uCar = car_vectors{ROI};
        vCar = car_vectors_2{ROI};
        
        %get smallest angle between vectors
        %theta(ROI) = atan2(uCar{1}(1)*uCar{2}(2) - uCar{1}(2)*uCar{2}(1), uCar{1}(1)*uCar{2}(1) + uCar{1}(2)*uCar{2}(2));
        theta(ROI) = atan2(uCar(1)*vCar(2) - uCar(2)*vCar(1), uCar(1)*vCar(1) + uCar(2)*vCar(2));
        theta_abs(ROI) = abs(theta(ROI));
    end
    
    %atan2(uCar(1)*vCar(2) - uCar(2)*vCar(1), uCar(1)*vCar(1) + uCar(2)*vCar(2));

    %% Verify a random ROI with a compass plot - looks right
    %plot unit vector of difference as well
if 0
    figure;
    for ii=1:100
        
        ROI = ii;
        %unit_vector_1 = exp(1i*ang_1(ROI)/(2*pi));
        %unit_vector_2 = exp(1i*ang_2(ROI)/(2*pi));
        diff_vector = exp(1i*theta_abs(ROI));
        
        compass(diff_vector,'k')
        hold on;
        compass(pf_vector_max{1}{4}(ROI),'r')
        compass(pf_vector_max{1}{5}(ROI),'r')
        pause
        clf
    end
end 
    
end


function [outputArg1,outputArg2] = centroid_diff_sessions(learn_comb_data,recall_comb_data,reg_learn, reg_recall)
%% Get all A/allB A&B


%% Set up variables

tuned_log_learning = learn_comb_data.tuned_log_learning;
pf_vector_max_learning = learn_comb_data.pf_vector_max_learning;


%% Make assignments based on tuning state 
%Calculate difference relative to D1 (for A, for B, for A vs.B - TS tuned/ min 5 events)

%get TS tuned - select only TS positive neurons in each reg column for A/B and A
%and B neurons

%LEARNING SETS
aa=2;
match_ROI.all = reg_learn{aa}.registered.multi.assigned_filtered(:,1:6);

%create separate A, B, AB matrices for TS filtering -duplicate
%match_ROI.all
match_ROI.allA = match_ROI.all;
match_ROI.allB = match_ROI.all;
match_ROI.AB = match_ROI.all;

%blank logical
blank_log = zeros(size(match_ROI.all,1),size(match_ROI.all,2));

for ss=1:6
    %get tuned indices for each session
    tuned_idx.allA = find(tuned_log_learning{aa}.tuned_logicals.tuned_log_filt_ts{ss}.allA  ==1 );
    [match_idx_A, column_match_idx{ss},~]  =  intersect(match_ROI.allA(:,ss),tuned_idx.allA);
    disp(length(match_idx_A))
    %translate match to logical value
    blank_log(column_match_idx{ss},ss) = 1;
    %nan non-matching ROIs
    match_ROI.allA(~blank_log(:,ss),ss) = nan;
end

%% Extract the tuning vectors for each class across days

%holder for each category vectors
max_vectors.A = nan(size(match_ROI.all,1),size(match_ROI.all,2));

%do for learning A first
for ss=1:6
    max_vectors.A(column_match_idx{ss},ss) = pf_vector_max_learning{aa}.pf_vector_max{ss}{4}( match_ROI.allA(column_match_idx{ss},ss));
end

%% Calculate anglular difference relative to d1 for all A

%get the vectors for both sessions
for ss=2:6
    compare_vectors.A{1,ss} = max_vectors.A(find(sum(~isnan(max_vectors.A(:,[1,ss])),2) ==2),[1,ss])
    
end

%convert complex vector to Cartesian coordinates
for ss=2:6
    for rr=1:size(compare_vectors.A{1, ss},1)
        %first ses
        compare_vec_cart.A{1,ss}{rr,1} = [real(compare_vectors.A{1, ss}(rr,1)), imag(compare_vectors.A{1, ss}(rr,1))];
        %second ses
        compare_vec_cart.A{1,ss}{rr,2} = [real(compare_vectors.A{1, ss}(rr,2)), imag(compare_vectors.A{1, ss}(rr,2))];
    end
end

for ss=2:6
    %calculate difference between each vector pair
    %for each ROI get difference between ts vector and reward vector
    for rr = 1:size(compare_vec_cart.A{1,ss},1)
        
        uCar = compare_vec_cart.A{1,ss}{rr,1};
        vCar = compare_vec_cart.A{1,ss}{rr,2};
        
        %get smallest angle between vectors
        theta{ss}(rr) = atan2(uCar(1)*vCar(2) - uCar(2)*vCar(1), uCar(1)*vCar(1) + uCar(2)*vCar(2));
        theta_abs{ss}(rr) = abs(theta{ss}(rr));
    end
end

%get mean and run paired wilcoxon
cellfun(@nanmean, theta_abs)


ranksum(theta_abs{2},theta_abs{3})

%plot cdf -test
figure
hold on
for ss=2:6
ecdf(theta_abs{ss})
pause
end


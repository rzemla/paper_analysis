function [recurr,frac_active] = recurrence_analysis(registered,options)


%set the number of sessions based on number of sessions imageed
nb_ses = size(options.sessionSelect,2);

%all active assigned
all_active_match = registered.multi.assigned_filtered;

%convert to binary
all_active_match_bin = ~isnan(all_active_match);

%SI binary match matrices for A and B trials
%SI
match_bin.si.A = ~isnan(registered.multi.matching_list_filtered.si_Aall_filt_event_filt);
match_bin.si.B = ~isnan(registered.multi.matching_list_filtered.si_Ball_filt_event_filt);
match_bin.si.AorB = match_bin.si.A | match_bin.si.B;

%TS
match_bin.ts.A = ~isnan(registered.multi.matching_list_filtered.ts_Aall_filt_event_filt);
match_bin.ts.B = ~isnan(registered.multi.matching_list_filtered.ts_Ball_filt_event_filt);
match_bin.ts.AorB = match_bin.ts.A | match_bin.ts.B;

%calculate recurrence
%for given 2 sessions, find neurons that are tuned in both divided by
%neurons that are tuned in either or both sessions
% (1/1)/(1/0 or 0/1 or 1/1)

%for session against the other
for ii=1:nb_ses
    for jj=1:nb_ses
        %TS A / B
        recurr.ts.A(ii,jj) =length(find(sum(match_bin.ts.A(:,[ii jj]),2) == 2))/...
            (length(find(sum(match_bin.ts.A(:,[ii jj]),2) == 1)) + length(find(sum(match_bin.ts.A(:,[ii jj]),2) == 2)));
        recurr.ts.B(ii,jj) =length(find(sum(match_bin.ts.B(:,[ii jj]),2) == 2))/...
            (length(find(sum(match_bin.ts.B(:,[ii jj]),2) == 1)) + length(find(sum(match_bin.ts.B(:,[ii jj]),2) == 2)));
        %TS A or B
        recurr.ts.AorB(ii,jj) =length(find(sum(match_bin.ts.AorB(:,[ii jj]),2) == 2))/...
            (length(find(sum(match_bin.ts.AorB(:,[ii jj]),2) == 1)) + length(find(sum(match_bin.ts.AorB(:,[ii jj]),2) == 2)));
        %SI A / B
        recurr.si.A(ii,jj) =length(find(sum(match_bin.si.A(:,[ii jj]),2) == 2))/...
            (length(find(sum(match_bin.si.A(:,[ii jj]),2) == 1)) + length(find(sum(match_bin.si.A(:,[ii jj]),2) == 2)));
        recurr.si.B(ii,jj) =length(find(sum(match_bin.si.B(:,[ii jj]),2) == 2))/...
            (length(find(sum(match_bin.si.B(:,[ii jj]),2) == 1)) + length(find(sum(match_bin.si.B(:,[ii jj]),2) == 2)));
        %SI A or B
        recurr.si.AorB(ii,jj) =length(find(sum(match_bin.si.AorB(:,[ii jj]),2) == 2))/...
            (length(find(sum(match_bin.si.AorB(:,[ii jj]),2) == 1)) + length(find(sum(match_bin.si.AorB(:,[ii jj]),2) == 2)));        
    end
end

%% Fraction of active neurons matched between sessions 
for ii=1:nb_ses
    for jj=1:nb_ses
        %relative to first session
        frac_active.first_ses(ii,jj) = length(find(sum(all_active_match_bin(:,[ii jj]),2) == 2))./length(find(sum(all_active_match_bin(:,ii),2) ==1));
        %relative to second session
        frac_active.second_ses(ii,jj) = length(find(sum(all_active_match_bin(:,[ii jj]),2) == 2))./length(find(sum(all_active_match_bin(:,jj),2) ==1));       
    end
end

%% Plot recurrence rate across time (relative to D1)
figure
hold on
title('T.S. tuned match rel D1')
plot(recurr.ts.A(1,:),'b')
plot(recurr.ts.B(1,:),'r')
%plot fraction active
ylim([0 1])
plot(frac_active.first_ses(1,:),'g')

figure
hold on
title('S.I. tuned match rel D1')
plot(recurr.si.A(1,:),'b')
plot(recurr.si.B(1,:),'r')
ylim([0 1])
plot(frac_active.first_ses(1,:),'g')

%A or B
figure
hold on
title('S.I. A or B tuned match rel D1')
plot(recurr.si.AorB(1,:),'m')
ylim([0 1])
plot(frac_active.first_ses(1,:),'g')

figure
hold on
title('T.S. A or B tuned match rel D1')
plot(recurr.ts.AorB(1,:),'m')
ylim([0 1])
plot(frac_active.first_ses(1,:),'g')

end


function [session_STCs] = multi_session_STC_extract(match_list_run_event_filtered,STC)
%INPUT (1): filtered list across days
%INPUT (2): chosen STC for each day 
%OUTPUT: cell with STCs matched on each session (nans where no match)

%nb of sessions
nb_ROI_list_match = size(match_list_run_event_filtered,1);

for ss=1:size(match_list_run_event_filtered,2)
    %get non-excluded positions for purposes of reconstruction
    reconstitute_idx{ss} = find(~isnan(match_list_run_event_filtered(:,ss)) ==1);
    %actual ROI indices to extract
    ROI_idx_temp = match_list_run_event_filtered(reconstitute_idx{ss},ss);
    
    %blank nan matrix - A trials
    session_STCs{1,ss} = nan(nb_ROI_list_match,100);
    
    %blank nan matrix - B trials
    session_STCs{2,ss} = nan(nb_ROI_list_match,100);
    
    %get ROI indices to extract on that session - A trials
    STC_sel{ss}.A = STC{ss}.A(:,ROI_idx_temp)';
    
     %get ROI indices to extract on that session - A trials
    STC_sel{ss}.B = STC{ss}.B(:,ROI_idx_temp)';   
    
    %insert STCs into correct match position - A
    session_STCs{1,ss}(reconstitute_idx{ss},:) = STC_sel{ss}.A;

    %insert STCs into correct match position - B
    session_STCs{2,ss}(reconstitute_idx{ss},:) = STC_sel{ss}.B;    
    
end


end


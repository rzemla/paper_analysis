function [traces_combined] = extract_traces(field_trace_events)
%extracts traces from each significant event to calculate mean


%combine traces into matrix to take mean
%length of each trace
trace_length = cellfun(@length,field_trace_events);
%max legnth of trace
max_length = max(cellfun(@length,field_trace_events));
%nb of traces
nb_traces= size(field_trace_events,2);
%matrix for combining traces
traces_combined = zeros(nb_traces,max_length);

%fill matrix with traces
for tt=1:nb_traces
    traces_combined(tt,1:trace_length(tt)) = field_trace_events{tt};
end



end


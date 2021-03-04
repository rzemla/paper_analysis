function [output_stats] = kruskalwallis_test(data_input)

%create numeric grouping vector for each class
group_vec = [];
for ii=1:numel(data_input)
    nb_neurons(ii) = numel(data_input{ii});
    group_vec = [group_vec; ii*ones(nb_neurons(ii),1)];
end

%insert data as column vector
[p, tbl, stats] = kruskalwallis(cell2mat(data_input),group_vec,'off');
%summary output data
output_stats = [p,tbl{2,5},stats.n,tbl{2,3}];

end


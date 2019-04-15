
function [on_map_shuffle_smooth]=gauss_filt(on_map_shuffle,Nbin,sigma)

for i=1:length(Nbin)
for n=1:size(on_map_shuffle{i},2)
on_map_shuffle_smooth{i}(:,n)=imgaussfilt(on_map_shuffle{i}(:,n),sigma);
end
end
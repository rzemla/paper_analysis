%sample 1-way data (LT TC TS A)
%test_data = [];

subject_data = repmat([1:3]', 1,5);

time_data = repmat([1:5],3,1);

%make matrix

comb_mat = [test_data(:),time_data(:), subject_data(:)];

%remove nans and split
%calculate GG correction (epsilon)
[x,y,gg] = GenCalcHFEps(test_data(:),[],[time_data(:)],subject_data(:))




%F dist params
df_n = 4;
df_d = 5;
f_value = 15.82;
%get p value
p_val = fcdf(f_value,df_n,df_d,'upper');

%with GG correction
epsilon_val = 0.001580;
%p_value with GG correction
p_val_gg = fcdf(f_value,df_n*epsilon_val,df_d*epsilon_val,'upper');


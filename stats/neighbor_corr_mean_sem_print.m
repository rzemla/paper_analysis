function [outputArg1,outputArg2] = neighbor_corr_mean_sem_print(ActXWord,WordHandle,txt_input, mean_both, sem_both, nb_comp)

%txt_input = 'PV Neighbor 1,2 vs. 6,7 A trial learning';
writeDefaultWordText(ActXWord,WordHandle,txt_input);
%newline
writeWordEnter(ActXWord,WordHandle,1);
%first mean and sem
mean_txt = mean_both(1);
sem_txt = sem_both(1);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

txt_input = ' vs. ';
writeDefaultWordText(ActXWord,WordHandle,txt_input);

%second mean and sem
mean_txt = mean_both(2);
sem_txt = sem_both(2);
write_mean_sem(ActXWord,WordHandle, mean_txt,sem_txt);

%newline
writeWordEnter(ActXWord,WordHandle,1);

%number of samples
txt_input = ['nb samp ', num2str(nb_comp)];
writeDefaultWordText(ActXWord,WordHandle,txt_input);
%newline
writeWordEnter(ActXWord,WordHandle,1);
writeWordEnter(ActXWord,WordHandle,1);

end


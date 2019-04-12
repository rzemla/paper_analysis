function [gaussFilter] = define_Gaussian_kernel(options)

%returns gaussian smoothing kernel of sigma defined by sigma filter

%%
%kernel size of Gaussian
gSigma = options.sigma_filter;
%half size of kernel
sz= 1 + 2*(3*gSigma);
%window vector
x_kernel = linspace(-sz / 2, sz / 2, sz);
%define filter as single term Gaussian
gaussFilter = exp(-x_kernel .^ 2 / (2 * gSigma ^ 2));
%normalize filter
gaussFilter = gaussFilter / sum (gaussFilter);

end


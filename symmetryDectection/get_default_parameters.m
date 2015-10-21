function parameters = get_default_parameters()
parameters = struct();
parameters.visualize = 0;
parameters.verbose = 1;

parameters.numberOfLines = 5;

parameters.searchRange = 2:2:80;
parameters.sigmas = 4;
parameters.searchAngles = [-pi/4, 0,pi/4];
parameters.symmetryAngles = ((1:16)/16)*pi;

parameters.filterCombinator = @(J1t, J2t) (J1t .* conj(J2t));
parameters.symmetryMetric = @(SIM) real(sqrt(SIM .* conj(SIM)));

parameters.miniumRhoDistance = 10;
parameters.segmentHistogramQuantileLow = .02;
parameters.segmentHistogramQuantileHigh = .92;

end
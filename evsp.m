%Quick little script to output solver result parameters
fprintf(' w_tot: %f \n Order NxM: %dx%d \n del: %f \n del t: %f \n EndTime: %.1f \n KxL elements: %dx%d \n Domain(WESN): %.3f %.3f %.3f %.3f \n TestCases: %d and %d more \n NearRange: %d \n Run time: %.0f \n',setup(1:13),length(setup)-15,setup(end-1:end))
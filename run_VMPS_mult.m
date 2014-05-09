for s = 1
	for delta = 0.05
		for alpha = 0.01:0.01:1
			try
				VMPS_FullSBM(s,alpha,delta)
            		catch err
                		getReport(err,'extended')
				fprintf('Problem with the combination: s=%.10g, alpha=%.10g, delta=%.10g',s,alpha,delta)
			end
		end
	end
end

exit

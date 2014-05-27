function run_VMPS_mult(p, eta)

for  period = p
    for factor = eta
		try
			VMPS_MLSBM(period,factor)
       	catch err
            getReport(err,'extended')
			fprintf('Problem with the combination: period=%.10g',period)
        end
    end
end

exit

end

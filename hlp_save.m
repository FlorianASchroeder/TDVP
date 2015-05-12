function hlp_save(fname,para,op,results,tresults,mps,Vmat)
% saves above variables serialized into a file.

para     = hlp_serialize(para);
save(fname,'para','-v7.3');
fprintf('.');

op       = hlp_serialize(op);
save(fname,'op','-append');
clear('op');
fprintf('.');

results  = hlp_serialize(results);
save(fname,'results','-append');
clear('results');
fprintf('.');

tresults = hlp_serialize(tresults);
save(fname,'tresults','-append');
clear('tresults');
fprintf('.');

mps      = hlp_serialize(mps);
save(fname,'mps','-append');
clear('mps');
fprintf('.');

Vmat     = hlp_serialize(Vmat);
save(fname,'Vmat','-append');
clear('Vmat');
fprintf('.\n');
end
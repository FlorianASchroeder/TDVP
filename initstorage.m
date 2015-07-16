function [op] = initstorage(mps, Vmat,op,para)
% op.Hlrstorage stores the sub hamiltonian of the (left?) and right block.
% op.Opstorage stores the operators in the (left?) and right bases.
% sweeps from right to left.
% Only use with rightnormalised MPS

% Commented by Florian Schroeder 29/01/2014
% Changed:
%   FS 20/10/2014: - use op = updateop() instead of explicit calculations
%                    to remove redundancy
%                  - moved separate treatment of L into for-loop

L 			  = length(mps);
M			  = size(op.h2term,1); 		% number of terms in sum per site
op.Hlrstorage = cell(1, L + 1);
op.Opstorage  = cell(M, 2, L+1);

% treat first and last terms seperately
op.Hlrstorage{1, 1} 	= 0;
op.Hlrstorage{1, L + 1} = 0;

% moved to updateop():
% op.Hlrstorage{L} = updateHright(op.Hlrstorage{L + 1}, op.h1term{L}, op.Opstorage(:,2,L+1), mps{L}, Vmat{L}, op.h2term(:,1,L), mps{L}, Vmat{L},M);
	% gives h1term{L} in effective basis of r_L-1 as:
	% op.Hlrstorage{l+1} = 0; op.h1term{L} = exists ; op.Opstorage(:,2,L+1) = [], op.h2term(:,1,L) = 0;
for m=1:M
	% moved to updateop():
%     op.Opstorage{m,2,L}=updateCright([],mps{L},Vmat{L},op.h2term{m,2,L},mps{L},Vmat{L});
	% transforms interaction terms into r_L-1 ef. basis
	op.Opstorage{m,1,1}	  =0;
	op.Opstorage{m,2,L+1} =0;
end

para.sweepto = 'l';

% middle terms build: l<-r
for j = L:-1:2
	para.sitej = j;
    op = updateop(op,mps,Vmat,j,para);
end

function [U, S, V, newDim] = truncateUSV(u,sv,v, para, minDim)
%% [U, S, V] = truncateUSV(U,S,V, para, minDim)
% Truncates Bond Dimensions according to Singular Values
% keeps always one SV in range [para.svmaxtol, para.svmintol]
% keeps always at least minDim dimensions
%
%	u : truncated in columns
%	s : vector of SV
%	v : is V', truncated in rows
%
% by Florian Schroeder 08/08/2015

if ~isempty(para)
	svmintol = para.svmintol;
	svmaxtol = para.svmaxtol;
else
	svmintol = 10^-4.5;
	svmaxtol = 10^1;
end

%% Truncate A dims
keepdims = find(sv > svmintol);
if length(keepdims) < minDim                     % keep at least Dmin bonds
	keepdims = 1:minDim;
end
% If smallest SV too large, keep 1 more
if (sv(keepdims(end)) > svmaxtol)
	if keepdims(end)+1 < length(sv)
		keepdims = [keepdims;keepdims(end)+1];
	else
		% expand Bond dims! Happens later?
		U = u; S = sv; V = v; newDim = length(sv);
		return;
	end
end
if length(keepdims) < length(sv)
	U = u(:,keepdims);                          % keep columns
	S = sv(keepdims);
	V = v(keepdims,:);                          % keep rows
    newDim = length(keepdims);
else
	U = u; S = sv; V = v; newDim = length(sv);
end

function conjGradient()
	%% how to optimise the truncated result?
end


end

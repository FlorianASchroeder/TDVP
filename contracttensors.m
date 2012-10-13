function [X, numindX] = contracttensors(X, numindX, indX, Y, numindY, indY)
%
Xsize = ones(1, numindX); Xsize(1:ndims(X)) = size(X);
Ysize = ones(1, numindY); Ysize(1:ndims(Y)) = size(Y);

indXl = 1:numindX; indXl(indX) = [];
indYr = 1:numindY; indYr(indY) = [];

sizeXl = Xsize(indXl);
sizeX = Xsize(indX);
sizeYr = Ysize(indYr);
sizeY = Ysize(indY);

if prod(sizeX)~= prod(sizeY)
    error('indX and indY are not of same dimension.');
end

if isempty(indYr)
    if isempty(indXl)
        X = permute(X, indX);
        X = reshape(X, [1, prod(sizeX)]);
        Y = permute(Y, indY);
        Y = reshape(Y, [prod(sizeY), 1]);
        X = X * Y; Xsize = 1;
        return;
    else
        X = permute(X, [indXl, indX]);
        X = reshape(X, [prod(sizeXl), prod(sizeX)]);
        Y = permute(Y, [indY]);
        Y = reshape(Y, [prod(sizeY), 1]);
        X = X * Y;
        Xsize = Xsize(indXl);
        X = reshape(X, [Xsize, 1]);
        return
    end
end

X = permute(X, [indXl, indX]); X = reshape(X, [prod(sizeXl), prod(sizeX)]);
Y = permute(Y, [indY, indYr]); Y = reshape(Y, [prod(sizeY), prod(sizeYr)]);

X = X * Y;
Xsize = [Xsize(indXl), Ysize(indYr)];
numindX = length(Xsize);
X = reshape(X, [Xsize, 1]);

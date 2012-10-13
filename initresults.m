function r=initresults(para)

L=para.L;
d_opt=para.d_opt;
D=para.D(end);

r.Vmat_vNE =zeros(1,L); %phonon basis von Neumann Amat_vNEj
r.Amat_vNE = zeros(1,L);
if para.parity~='n'
    r.Vmat_sv=cell(1,L);
    %r.Vmat_sv=zeros(L,max(d_opt)/2);
    r.Amat_sv=cell(1,L-1);
    %r.Amat_sv=zeros(L,D/2);
else
    r.Vmat_sv=cell(1,L);
    %r.Vmat_sv=zeros(L,max(d_opt));
    r.Amat_sv=cell(1,L-1);
    %r.Amat_sv=zeros(L,D);
end

r.leftge = zeros(L,1);
r.geoffset=zeros(L,1);
r.Eerror(1)=1;
r.lastVmat_vNE=zeros(1,L);
end
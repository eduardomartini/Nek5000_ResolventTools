function G2 = interpTensors(f1,G1,f2,defValue)
%  G2 = interpTensors(f1,G1,f2,defValue)
% Interpolates a tensor G(i,j,it), with frequencies given by f1 in the last
% entry onto frequencies f2
% If a frequency is not in the range of f2, defValue is used 
if ~exist('defValue'); defValue = 0; end
G2 = zeros(size(G1,1),size(G1,2),length(f2));
for i=1:size(G1,1)
    for j = 1:size(G1,2)
        yinterp(:) = G1(i,j,:);
        G2(i,j,:)= fftshift(interp1(ifftshift(f1),ifftshift(yinterp),ifftshift(f2)));
    end
end

for i=1:length(f2)
    if sum(isnan(G2(:,:,i)))>0
        G2(:,:,i) = defValue;
    end
end
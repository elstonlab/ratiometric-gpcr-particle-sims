function [COMvec,t,nTarg,confidence,normProj,xyzAlignment,xyAlignment]=quantify_v2(filename,nframes,targfield,totalTarg)
% We only use the xyAlignment (the angle between the +x direction and the center of mass vector
% for later analysis. However, many other quantities (such as the ones provided here) can be
% computed.
%
% Input: filename: source file from smoldyn output (.xyz) to quantify.
%        nframes: expected # frames in dataset (typically 62)
%        targfield: 'name' of target field. typically 'Ga' to quantify active G protein
%        totalTarg: total amount of molecule, inactive+active. Typically 2500 for G protein.
%
% Outputs:
%        COMvec: center of mass vector
%        t: time
%        nTarg: # of target protein
%        confidence: confidence (Lakhani and Elston, 2017)
%        normProj: normalized projection (dot product normalized to magnitude)
%        xyzAlignment: xyz alignment with gradient, arccos of normalized projection (degrees)
%        xyAlignment: xy alignment with gradient after discarding z; arccos of normalized [xy only] projection (deg)
%
% Aside: nframes is typically 62 or 602 frames. However, the data generally had
% only 61 or 601 "real" frames (0,10,...,600) or (0,1,2,...,600). 
% Due to an old typo in my configuration files, the t0 is duplicated.
% As a result, we ignore the first frame in later analyses.

[t,positions]=read_molPos3(filename,nframes);
COMvec=zeros(nframes,3);
nTarg=zeros(nframes,1);
confidence = zeros(nframes,1);
normProj = zeros(nframes,1);
xyzAlignment = zeros(nframes,1);
xyAlignment = zeros(nframes,1);

for i=1:nframes
   if isnan(t(i))
       break;
   end
   nTarg(i) = numel(positions.(targfield){i}(:,1));
   COMvec(i,:)=mean(positions.(targfield){i},1);
   
   xyAlignment(i) = rad2deg(acos(normalized_projection(COMvec(i,1:2)-2.5,[1 0])));
   
   confidence(i) = 1/2.5.*(nTarg(i)/totalTarg).*sqrt(sum((COMvec(i,:)-2.5).^2));
end
end

function norm_proj = normalized_projection(v1,v2)
% Input: vectors 1 and 2 of same length.
% Output: the normalized projection of vector 1 onto vector 2 (or vice-versa)
%         0 means no colinearity at all.
%        -1 is perfect anti-colinearity
%        +1 is perfect colinearity
norm_proj = dot(v1,v2)/(norm(v1)*norm(v2));
end
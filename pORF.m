function prob=pORF(N_seq,N_orf)

ORFproteinLength = ceil(N_orf/3);
for frame = 1:3 %Loop over the frame
proteinLength = floor((N_seq-frame+1)/3); %get length of protein sequence
allprobs = zeros(proteinLength,1);
for pos = 1:proteinLength
    allprobs(pos) = pORFLongerThanOmin(pos,ORFproteinLength,proteinLength);
end
pNoORF(frame) = prod(-1*allprobs+1);
end
prob = 1-prod(pNoORF);




function prob=pORFLongerThanOmin(pos,Omin,N_seq)
ORFlengths = Omin:(N_seq-pos+1); %possible lengths of ORF
if isempty(ORFlengths)
    prob = 0;
else
    allprobs = arrayfun(@pORFLength,ORFlengths);
    prob = sum(allprobs);
end


function prob = pORFLength(k)
%probability that a sequence of length k is an ORF.
%Must begin with start, end with stop and all middle
%codons must be non-stop.
%note: k is length of protein sequence. 
Pstart = 1/64;
Pstop = 3/64;
PNotStop = 1-Pstop;
prob=Pstart*Pstop*PNotStop^(k-2);
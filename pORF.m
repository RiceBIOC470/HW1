function prob=pORF(N_seq,N_orf)
%Computes the probability that there an ORF or length N_orf or longer in
%DNA Sequence of length DNA.

%This solution is close to exact but makes one simplification: it does
%not consider correlations between the reading frames. For example, if
%there is a start codon at position 1, there cannot also be one a position
%two or three. Conversely if there no start codon at 1, it is slightly more
%likely to have one at two or three. Considering these correlations
%increases the difficulty of the problem considerably without much added
%accuracy (as you can see from the plot in the solutions). 

ORFproteinLength = ceil(N_orf/3); %convert ORF length to a length of protein sequence
for frame = 1:3 %Loop over the frame each is a different protein sequence
    proteinLength = floor((N_seq-frame+1)/3); %get length of protein sequence
    allprobs = zeros(proteinLength,1);
    for pos = 1:proteinLength
        %this computes the probability than an ORF of at least the right length
        %begins at each position.
        allprobs(pos) = pORFLongerThanOmin(pos,ORFproteinLength,proteinLength);
    end
    %1-product(allprobs) gives the probability that there are no ORFs long enough in this
    %frame
    pNoORF(frame) = prod(-1*allprobs+1);
end
%prod(pNoORF) gives the probabilty that there is no ORF long enough in any
%frame. 1-prod(pNoORF) gives the probability that there is one.
prob = 1-prod(pNoORF);

function prob=pORFLongerThanOmin(pos,Omin,N_seq)
%This is the probability that in a sequence of length N_seq
%an ORF of length greater than or equal to Omin begins at position pos.
%Note: these are amino acid sequence lengths and positions, not DNA bases.
ORFlengths = Omin:(N_seq-pos+1); %possible lengths of ORF
if isempty(ORFlengths)
    prob = 0;
else
    %the function arrayfun
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
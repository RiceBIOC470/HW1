% HW1 challenge problems

% 1. Write a solution to Problem 2 part 2 that doesn't use any loops at
% all. Hint: start by using the built in function bsxfun to make a matrix of all distances
% between start and stop codons. 


N = 500;

dnaletters = ['A' 'T' 'G' 'C'];
rand_seq = dnaletters(randi(4,1,N));

startcodon_pos = strfind(rand_seq, 'ATG');
stopcodon_pos = [strfind(rand_seq, 'TAA') strfind(rand_seq, 'TGA') strfind(rand_seq, 'TAG')];

found = 0; % variable used to check if any ORF is found in the sequence. 

if (~isempty(startcodon_pos) && ~isempty(stopcodon_pos))
    orflength = bsxfun(@minus, stopcodon_pos, startcodon_pos');
    % orflength is a matrix with lengths of all dna sequences starting with start codon and ending with stop codon.
    % number of rows  in orflength is the same as number of start_codons, number of columns equals the number of stop codons.
    
    condition1 = orflength > 0; % ensures stop codon is after start codon
    condition2 = mod(orflength,3) == 0; % ensures that the start and stop codons are in frame
    
    orflength(~(condition1&condition2)) = 0; % assigns zero values for  orfs that do not satisfy the above two conditions
    
    if(sum(sum(orflength)) > 0) % checks if any orf is found
        orflength(orflength==0) = NaN; % filters out zero values.
        [allorfs_length, allorfs_stop] = min(orflength, [], 2); % finds the length and stop codon position of the first in frame orf for each start codon.
        orfstart = find(allorfs_length == max(allorfs_length)); % finds start codon corresponding to longest in frame orf.
        orfstop = allorfs_stop(orfstart); % gives stop codon corresponding to longest in frame orf.
        found = 1;
        
             
    end
end


if(found ==0)
    display('No ORF found');
else
    orf_sequence = rand_seq(startcodon_pos(orfstart): stopcodon_pos(orfstop)+2); % gives the sequence of the orf.
    orf_length = length(orf_sequence);
    
    
    display(['Longest ORF found:', orf_sequence]);
    display(['Length:', int2str(orf_length)]);
       
end

%%

% 2. Problem 2, part 4. Use Matlab to compute the exact solution to this
% problem and compare your answer to what you got previously by testing
% many sequences. Plot both on the same set of axes. Hint: to get started 
% think about the following:
% A. How many sequences of length N are there?
% B. How many ways of making an ORF of length N_ORF are there?
% C. For each N_ORF how many ways of position this reading frame in a
% sequence of length N are there?

% 3. Problem 3. Assume that the error in each Cp is the standard deviation
% of the three measurements. Add a section to your code that propogates this
% uncertainty to the final results. Add error bars to your plot. (on
% propagation of error, see, for example:
% https://en.wikipedia.org/wiki/Propagation_of_uncertainty

genestoanalyse = 1:3;
fold_change = zeros(3,5); % A matrix to store fold change of genes in 5 conditions.
standard_deviation = zeros(3,5);% A matrix to store the std deviation of fold change values for three genes.

gene_number = 1;

for i = 1:3:9 %column numbers corresponding to different genes in the variable platedata.
    for j = 2:6 %row numbers corersponding to different conditions in the variable platedata.
        
        Cp0 = mean(platedata(1,i:i+2)); % platedata is a matrix containing Cp values in the desired format Problem 3 part 2. 
        CpX = mean(platedata(j,i:i+2)); 
        
        CpN0 = mean(platedata(1, 10:12));
        CpNX = mean(platedata(j, 10:12));
        
        fold_change(gene_number, j-1) = 2^[Cp0 - CpX - (CpN0 - CpNX)];
        
        for k = 1:3 %computes the fold change for each gene triplicate to find the standard deviation across them. 
            cp0 = platedata(1,i+(k-1)); % The Cp value of one of the triplicates(i + (k-1)) of gene i in condition 1
            cpx = platedata(j,i+(k-1)); % The Cp value that triplicate in condition j
            
            foldchange_1(k) = 2^[cp0 - cpx - (CpN0 - CpNX)];
        end
        
        standard_deviation(gene_number, j-1) = std(foldchange_1); % The function "std" computes the standard deviation for across the three measurements stored in foldchange_1.
        
    end
    gene_number = gene_number + 1;
end



figure;
h = bar(fold_change);
hold on;

for i = 1:5
    hold on;
    errorbar((1:3) + h(i).XOffset, (fold_change(:,i))', (standard_deviation(:,i))', 'k.');
end

ylabel('Fold Change');

label = [{'Condition 1'}, {'Condition 2'}, {'Condition 3'}, {'Condition 4'}, {'Condition 5'}];
legend(label);

ax = gca;
ax.XTickLabel = [{'Gene 1'}, {'Gene 2'}, {'Gene 3'}];
ax.FontSize = 16;


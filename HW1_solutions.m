% Homework 1. Due before class on 9/6/16.

%% Problem 1 - addition with strings

% Fill in the blank space in this section with code that will add
% the two numbers regardless of variable type. Hint see the matlab
% functions ischar, isnumeric, and str2num.

%your code should work no matter which of these lines is uncommented.
%x = 3; y = 5; % integers
x = '3'; y= '5';

%your code goes here
if ischar(x)  %if it is character, convert it to number
    x = str2num(x);
end

if ischar(y) %same for y
    y = str2num(y);
end

sumvalue = x+y;

%output your answer
display(['x + y = ' int2str(sumvalue)]);
%% Problem 2 - our first real biology problem. Open reading frames and nested loops.

%part 1: write a piece of code that creates a random DNA sequence of length
% N (i.e. consisting of the letters ATGC) where we will start with N=500 base pairs (b.p.).
% store the output in a variable
% called rand_seq. Hint: the function randi may be useful. Bonus points if
% you can do this without using a loop.

N = 500; % define sequence length
dnaletters = ['A', 'T', 'G', 'C']; %define the letters
rand_seq = dnaletters(randi(4, 1, N));
% randi generates a random sequence of numbers having the values 1 to 4.
% These values are then used to refer to the 4 characters (A,T,G,C) of the
% string dnaletters to generate a DNA sequence

%part 2: open reading frames (ORFs) are pieces of DNA that can be
% transcribed and translated. they start with a start codon (ATG) and end with a
% stop codon (TAA, TGA, or TAG). Write a piece of code that finds the longest ORF
% in your seqeunce rand_seq. Hint: see the function strfind.

startcodon_pos = strfind(rand_seq, 'ATG'); % starting position of start codon in the sequence.
stopcodon_pos = [strfind(rand_seq, 'TAA') strfind(rand_seq, 'TAG') strfind(rand_seq, 'TGA')]; % starting position of stop codon in the sequence

if ~isempty(startcodon_pos) && ~isempty(stopcodon_pos)
    
    firstStopCodon = zeros(1,length(startcodon_pos)); %start will all zeros
    for ii = 1:length(startcodon_pos) % loop over start codons
        ORFlengths = stopcodon_pos - startcodon_pos(ii); %distance from current start codons to all stop codons
        good_length = 1e8; %store the length here. initialize it to something high, so can only go down when ORF is found
        good_index = 0; % will store the index of the ORF here
        for jj = 1:length(ORFlengths)
            if ORFlengths(jj) > 0 && ... %stop after start
                    mod(ORFlengths(jj),3) == 0 && ... %its divisible by 3
                    ORFlengths(jj) < good_length %its shorter than the current one found
                good_length = ORFlengths(jj);
                good_index = jj;
            end
        end
        if good_index > 0 % check if something was found
            firstStopCodon(ii) = stopcodon_pos(good_index); %if so, store it
        else
            firstStopCodon(ii) = startcodon_pos(ii); %ORF end is the same as its beginning, length is zero
        end
    end
    ORFSizes = firstStopCodon - startcodon_pos + 3; %get the sizes of all ORFS
    [ORFlengthMax, ind_max] = max(ORFSizes); % take the maximum
end

if exist('ORFlengthMax', 'var') &&  ORFlengthMax > 0 %if something was found display the results.
    disp(['Longest open reading frame is ' int2str(ORFlengthMax) ' base pairs long. Start codon at '...
        int2str(startcodon_pos(ind_max)) '. Stop codon at ' int2str(firstStopCodon(ind_max))]);
else
    disp('No ORF found');
end
%%
%part 3: copy your code in parts 1 and 2 but place it inside a loop that
% runs 10000 times. Use this repeating operation to determine the probability
% that an sequence of length 500 has an ORF of greater than 50 b.p.

counter = 0; % counter increases by 1 if an ORF is found.

for tt = 1:10000
    
    N = 500; % define sequence length
    dnaletters = ['A', 'T', 'G', 'C'];
    rand_seq = dnaletters(randi(4, 1, N));
    
    startcodon_pos = strfind(rand_seq, 'ATG'); % starting position of start codon in the sequence.
    stopcodon_pos = [strfind(rand_seq, 'TAA') strfind(rand_seq, 'TAG') strfind(rand_seq, 'TGA')]; % starting position of stop codon in the sequence
    
    if ~isempty(startcodon_pos) && ~isempty(stopcodon_pos)
        
        firstStopCodon = zeros(1,length(startcodon_pos)); %start will all zeros
        for ii = 1:length(startcodon_pos) % loop over start codons
            ORFlengths = stopcodon_pos - startcodon_pos(ii); %distance from current start codons to all stop codons
            good_length = 1e8; %store the length here. initialize it to something high, so can only go down when ORF is found
            good_index = 0; % will store the index of the ORF here
            for jj = 1:length(ORFlengths)
                if ORFlengths(jj) > 0 && ... %stop after start
                        mod(ORFlengths(jj),3) == 0 && ... %its divisible by 3
                        ORFlengths(jj) < good_length %its shorter than the current one found
                    good_length = ORFlengths(jj);
                    good_index = jj;
                    
                end
            end
            if good_index > 0 % check if something was found
                firstStopCodon(ii) = stopcodon_pos(good_index); %if so, store it
                
            else
                firstStopCodon(ii) = startcodon_pos(ii); %ORF end is the same as its beginning, length is zero
            end
        end
        ORFSizes = firstStopCodon - startcodon_pos + 3; %get the sizes of all ORFS
        
        if sum(ORFSizes > 3) > 0  && sum(ORFSizes > 50) > 0 % checks if an ORF longer than 50 is found
            counter = counter +1;
        end
        
    end
end

prob = counter/1e4;
display(['Probability of finding an ORF longer than 50bp: ' num2str(prob)]);
%%
%part 4: copy your code from part 3 but put it inside yet another loop,
% this time over the sequence length N. Plot the probability of having an
% ORF > 50 b.p. as a funciton of the sequence length.
sequence_lengths = 0:50:1000;
nseq = length(sequence_lengths);

prob = zeros(nseq,1);
for kk = 1:nseq
    
    N = sequence_lengths(kk);

    counter = 0; % counter increases by 1 if an ORF is found.
    
    for tt = 1:10000
        
        % define sequence length
        dnaletters = ['A', 'T', 'G', 'C'];
        rand_seq = dnaletters(randi(4, 1, N));
        
        startcodon_pos = strfind(rand_seq, 'ATG'); % starting position of start codon in the sequence.
        stopcodon_pos = [strfind(rand_seq, 'TAA') strfind(rand_seq, 'TAG') strfind(rand_seq, 'TGA')]; % starting position of stop codon in the sequence
        
        if ~isempty(startcodon_pos) && ~isempty(stopcodon_pos)
            
            firstStopCodon = zeros(1,length(startcodon_pos)); %start will all zeros
            for ii = 1:length(startcodon_pos) % loop over start codons
                ORFlengths = stopcodon_pos - startcodon_pos(ii); %distance from current start codons to all stop codons
                good_length = 1e8; %store the length here. initialize it to something high, so can only go down when ORF is found
                good_index = 0; % will store the index of the ORF here
                for jj = 1:length(ORFlengths)
                    if ORFlengths(jj) > 0 && ... %stop after start
                            mod(ORFlengths(jj),3) == 0 && ... %its divisible by 3
                            ORFlengths(jj) < good_length %its shorter than the current one found
                        good_length = ORFlengths(jj);
                        good_index = jj;
                        
                    end
                end
                if good_index > 0 % check if something was found
                    firstStopCodon(ii) = stopcodon_pos(good_index); %if so, store it
                    
                else
                    firstStopCodon(ii) = startcodon_pos(ii); %ORF end is the same as its beginning, length is zero
                end
            end
            ORFSizes = firstStopCodon - startcodon_pos + 3; %get the sizes of all ORFS
            
            if sum(ORFSizes > 3) > 0  && sum(ORFSizes > 50) > 0 % checks if an ORF longer than 50 is found
                counter = counter +1;
            end
            
        end
    end
    
    
    prob(kk) =  counter/1e4 ;
end

figure;
plot(sequence_lengths, prob);
xlabel('Sequence Length');
ylabel('Probability');

ax = gca; % gets the axis information of the current figure.
ax.FontSize = 14; % sets the FontSize of the text in the figure to the assigned value.

%part 5: Make sure your results from part 4 are sensible. What features
% must this curve have (hint: what should be the value when is small or when
% N is very large? how should the curve change in between.) Make sure your
% plot looks like this.
% --- Ans) The probability of finding an ORF greater than 50bp is 0 for sequences less than 50bp.
% It increases as the value of N increases, and finally saturates to 1.

%% problem 3 data input/output and simple analysis

%The file qPCRdata.txt is an actual file that comes from a Roche
%LightCycler qPCR machine. The important columns are the Cp which tells
%you the cycle of amplification and the position which tells you the well
%from the 96 well plate. Each column of the plate has a different gene and
%each row has a different condition. Each gene in done in triplicates so
%columns 1-3 are the same gene, columns 4-6 the same, etc.
%so A1-A3 are gene 1 condition 1, B1-B3 gene 1 condition 2, A4-A6 gene 2
%condition 1, B4-B6 gene2 condition 2 etc.

% part1: write code to read the Cp data from this file into a vector. You can ignore the last two
% rows with positions beginning with G and H as there were no samples here.

file1 = fopen('qPCRdata.txt', 'r'); % for this command to work, make sure the file is in the current working directory. Otherwise, specify the full path of the file.

cp_column = 5; % column number corresponding to Cp values.
pos_column = 3; % column number corresponding to position values.
firstrow_found = 0; % variable used to check if line corresponding to A1 position is found.

cp_values = 0; % variable that stores all the cp_values for positions A1 to F12.
counter = 1; % variable that increases by 1 for every new Cp value added to the array cp_values.


while file1 ~= -1
    row1 = fgets(file1); %gets a new line of file file1 with every iteration of the loop. Stores information as a string.
    row1_split = strsplit(row1, '\t'); % splits the string into a substrings whenever a tab is encountered. Each substring represents a different column in the line.
    
    
    if numel(row1_split) >= cp_column  % checks if the row has at least as many columns that can include cp values
        position = row1_split{pos_column}; % gets the value of position.
        
        if firstrow_found == 0   % checks if the row corresponding to A1 is not found yet.
            if strcmp(position, 'A1') % checks if the position is A1.
                firstrow_found = 1; % modifies the value of firstrow_found.
                cp_values(counter) = str2num(row1_split{cp_column}); %converts the string in column corresponding to Cp values into integer and assigns it to the first value of array cp_values.
                counter = counter+1; %increments counter.
            end
        else
            if strcmp(position, 'G1') % breaks the while loop i.e. ends the loop if we have reached the line corresponding to position 'G1'.
                break;
            else
                cp_values(counter) = str2num(row1_split{cp_column});
                counter = counter+1;
            end
            
            
        end
    end
end

%%
% Part 2: transform this vector into an array representing the layout of
% the plate. e.g. a 6 row, 12 column array should that data(1,1) = Cp from
% A1, data(1,2) = Cp from A2, data(2,1) = Cp from B1 etc.

platedata = zeros(6,12); %the matrix platedata refers to the array that represents the layout of the plate.
row_num = 1; % denotes the rows of the matrix platedata, starting with row1.

% In the loop that follows, i values correspond to the position values (A1-A12....F1-F12). In the matrix platedata, each row contains information for 1 condition.
% Conditions change every 12 positions. With every iteration of the loop, a different condition is called and the values are assigned to a new row in the matrix platedata.
for i = 1:12:72
    platedata(row_num,1:12) = cp_values(i:i+11);
    row_num = row_num+1;
end

%%
% Part 3. The 4th gene in columns 10 - 12 is known as a normalization gene.
% That is, it should not change between conditions and it is used to normalize
% the expression values for the others. For the other three
% genes, compute their normalized expression in all  conditions, normalized to condition 1.
% In other words, the fold change between these conditions and condition 1. The
% formula for this is 2^[Cp0 - CpX - (CpN0 - CpNX)] where Cp0 is the Cp for
% the gene in the 1st condition, CpX is the value of Cp in condition and
% CpN0 and CpNX are the same quantitites for the normalization gene.
% Plot this data in an appropriate way.

genestoanalyse = 1:3;
fold_change = zeros(3,5); % A matrix to store fold change of genes in 5 conditions. Each row represents a gene and each column represents a condition.

gene_number = 1; % variable used to refer to the rows of matrix fold_change.

for i = 1:3:9 %column numbers corresponding to different genes in the matrix platedata. Genes change every 3 columns.
    Cp0 = mean(platedata(1,i:i+2)); % mean across gene triplicates to get the Cp value of gene in column i in condition 1
    CpN0 = mean(platedata(1, 10:12)); % mean across gene triplicates to get the Cp value of normalization gene in condition 1    
    for j = 2:6 %row numbers corersponding to different conditions in the matrix platedata.
        
        CpX = mean(platedata(j,i:i+2)); % mean across gene triplicates to get the Cp value of gene in column i in condition j
        
        CpNX = mean(platedata(j, 10:12)); % mean across gene triplicates to get the Cp value of normalization gene in condition j
        
        fold_change(gene_number, j-1) = 2^[Cp0 - CpX - (CpN0 - CpNX)]; % fold change for each gene stored in matrix fold_change.
        
        
    end
    gene_number = gene_number + 1; %increments the gene_number.
end

%
figure;
bar(fold_change); %plots a grouped bar graph to represent fold_change. It is plotted as three groups of 5 bars each(corresponding to fold_change in five different conditions). Each group corresponds to a different gene.
ylabel('Fold Change');

label = [{'Condition 1'}, {'Condition 2'}, {'Condition 3'}, {'Condition 4'}, {'Condition 5'}];
legend(label);

ax = gca;
ax.XTickLabel = [{'Gene 1'}, {'Gene 2'}, {'Gene 3'}];
ax.FontSize = 16;


%% problem 4. starting with image data:

%download this zstack of a fly embryo with all the cell membranes marked:
%https://www.dropbox.com/s/atq9j35gtq6jcov/Screenshot%202016-08-25%2013.40.53.png?dl=0

% part 1. each image can be read in with the imread command e.g. img =
% imread('z1.tif'); will store the data in img for the first image.
% wrote a loop that reads all the data and makes a maximum intensity
% project (that is, each pixel in the image should contain the maximal
% value for that pixel in the entire z-stack. store this in the variable
% max_img


totalimages = 384; % number of '.tif' files in the folder flydata.

for ii = 1:totalimages %loop over all files
    image_name = ['flyData/z' int2str(ii) '.tif']; % name of the file.
    %You can also use an alternate command: image_name = sprintf('z%d.tif', ii). 
    %The function sprintf replaces the value of ii in the place of %d in the string 'z%d.tif' . 
    %If ii = 1, the output will be z1.tif.
    
    im1 = imread(image_name); %reads the image file.
    
    if(ii==1)
        max_img = im1; %assigns max_img to image1 for the first image
    else
        %for each subsequent image, compares the value of every pixel in 
        %max_img with the new image read (im1) and keeps the maximum 
        %value of the two.
        max_img = max(max_img, im1); 
    end
end

% part 2. write a loop that reads all the images and writes them as a
% single 3D image. hint: see the function imwrite and its option 'append'

new3dfile_name = 'max_img_3d.tif'; % name of the new file.

for ii = 1:totalimages
    image_name = ['flyData/z' int2str(ii) '.tif'];
    im1 = imread(image_name);
    
    if(ii==1)
        imwrite(im1, new3dfile_name); % writes the first image to the new file.
    else
        imwrite(im1, new3dfile_name, 'WriteMode', 'append'); % appends every subsequent image to the new file
    end
end
%%
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

% Solution: the code pORF.m in this folder contains the function
% pORF(N_Seq,N_ORF) computes the likelihood that a sequence of legnth N_Seq
% contains an ORF of length at least N_ORF. This calls that function for a
% number of different sequence lengths:
sequence_lengths = 0:50:1000;
nseq = length(sequence_lengths);
prExact = zeros(nseq,1);
for jj=1:nseq
    prExact(jj) =pORF(sequence_lengths(jj),50);
end

%above the approximate solution by testing repeated sequences is stored in 
%the variable prob. This compares the two. 

% make a plot to compare the 2:
figure; hold on;
plot(sequence_lengths,prExact,'k-','LineWidth',3);
plot(sequence_lengths,prob,'rs','MarkerSize',18);
xlabel('Sequence Length','FontSize',24);
ylabel('Probability','FontSize',24);
legend({'Exact Solution','Testing sequences'},'FontSize',18);

%%
% 3. Problem 3. Assume that the error in each Cp is the standard deviation
% of the three measurements. Add a section to your code that propogates this
% uncertainty to the final results. Add error bars to your plot. (on
% propagation of error, see, for example:
% https://en.wikipedia.org/wiki/Propagation_of_uncertainty

fold_change_error = zeros(3,5); % A matrix to store the error in the fold change of genes in 5 conditions.
gene_number = 1;

% calculating fold_change error.

% std(F)/F = log(2)*sqrt[(std(Cp0)/Cp0)^2 + (std(CpX)/CpX)^2) + (std(CpN0)/CpN0)^2 + (std(CpNX)/CpNX)^2].
% Error in fold change is represented as its std deviation - std(F).
% Since fold change in each gene in each condition is the function of Cp
% values of that gene and the normalisation gene in condition 1, the error
% in fold change will take into account the standard deviation in each of
% these four values. 
data = platedata;
for i = 1:3:9 %column numbers corresponding to different genes in the variable data.
    for j = 2:6 % 
        
        Cp0 = mean(data(1,i:i+2)); % taking mean across gene triplicates to get the Cp value of gene i in condition
        std0 = std(data(1,i:i+2)); % taking std deviation across gene triplicates. 
        
        CpX = mean(data(j,i:i+2));
        stdX = std(data(j,i:i+2));
        
        CpN0 = mean(data(1, 10:12));
        stdN0 = std(data(1,10:12));
        
        CpNX = mean(data(j, 10:12));
        stdNX = std(data(j,10:12));
        
        all_reqvalues = [std0 stdX stdN0 stdNX]; %this variable defined just to make calculation in next step less clumsy.
 
        fold_change_error(gene_number, j-1) = log(2)*fold_change(gene_number, j-1)*sqrt(sum(all_reqvalues.^2));
    end
    
    gene_number = gene_number + 1;
end
        
% Make the plot
figure;
h = bar(fold_change);
hold on;

%errorbar makes the plot. Notice the use of the XOffset property in the bar
%plot to get the errorbars to line up correctly. 
for i = 1:5
    hold on;
    errorbar((1:3) + h(i).XOffset, (fold_change(:,i))', (fold_change_error(:,i))', 'k.');
end

ylabel('Fold Change');

label = [{'Condition 1'}, {'Condition 2'}, {'Condition 3'}, {'Condition 4'}, {'Condition 5'}];
legend(label);

ax = gca;
ax.XTickLabel = [{'Gene 1'}, {'Gene 2'}, {'Gene 3'}];
ax.FontSize = 16;




library(VLF) # load VLF package

### VLF analysis on full birds dataset, N = 11333 barcodes ###

data(birds) # load birds nucleotide dataset

out_birds <- vlfFun(birds)

sum(out_birds$specimen) # total number of specimen VLFs - 771
sum(out_birds$position) # total number of positionsl VLFs - 771

length(which(out_birds$specimen != 0)) # total number of specimens with VLFs - 552
length(which(out_birds$position != 0)) # # total number of positions with VLFs - 241

which(out_birds$specimen != 0) # specimens with VLFs
which(out_birds$position != 0) # sites containing VLFs

rowSums(out_birds$sas) # number of singleton and shared VLFs

species.names <- birds[,2]
specimen.Number <- nrow(birds)
rownames(birds) <- species.names
Nuc.count <- count.function(birds, specimen.Number, seqlength = 648)
frequency.matrix <- ffrequency.matrix.function(Nuc.count, seqlength = 648)
birdSpec.freq <- specimen.frequencies(frequency.matrix, birds, specimen.Number, species.names, seqlength = 648)
Bird_specimen_VLFcount <- VLF.count.spec(birdSpec.freq, p = 0.001, seqlength = 648)
Bird_position_VLFcount <- VLF.count.pos(birdSpec.freq, p = 0.001, seqlength = 648)
bird_VLFconvert <- VLF.convert.matrix(birds, birdSpec.freq, p = 0.001, seqlength = 648)
bird_VLFnuc <- VLF.nucleotides(bird_VLFconvert, birds, seqlength = 648)
bird_VLFreduced <- VLF.reduced(bird_VLFnuc, Bird_specimen_VLFcount, seqlength = 648)
bird_species <- separate(bird_VLFreduced)
birds_singleAndShared <- find.singles(bird_species, seqlength = 648)

# singleton error rate
birds_single_error <- Error.Rate(birds_singleAndShared[1,], spec = specimen.Number, seqlength = 648) 
# shared error rate
birds_shared_error <- Error.Rate(birds_singleAndShared[2,], spec = specimen.Number, seqlength = 648)
# total error rate
birds_error_total <- Error.Rate(birds_singleAndShared[1,] + birds_singleAndShared[2,], spec = specimen.Number, seqlength = 648)

Decile.Plot(birds_singleAndShared, seqlength = 648) # plot decile plot
Sliding.Window(birds_singleAndShared, seqlength = 648, n = 30) # plot sliding window


### VLF analysis for Canada goose (Branta canadensis), N =  125 barcodes ###

x <- separate(birds) # partition into lists by species name

out_goose <- vlfFun(birds, own = x[[317]]) # B. canadensis is list element 317

rowSums(out_goose$sas) # number of singleton and shared VLFs


sum(out_goose$ownSpecCount) # total number of specimen VLFs - 27
sum(out_goose$ownPosCount) # total number of positionsl VLFs - 27

length(which(out_goose$ownSpecCount != 0)) # total number of specimens with VLFs - 18
length(which(out_goose$ownPosCount != 0)) # total number of positions with VLFs - 10

which(out_goose$ownSpecCount != 0) # specimens with VLFs 
which(out_goose$ownPosCount != 0) # sites containing VLFs

table(out_goose$ownSpecCount) # frequency of specimen VLFs
table(out_goose$ownPosCount) # frequency of positional VLFs

goose_single_error <- Error.Rate(out_goose$sas[1,], out_goose$sas[2,], spec = 125, seqlength = 648) 
goose_shared_error <- Error.Rate(out_goose$sas[2,], out_goose$sas[1,], spec = 125, seqlength = 648)
goose_totalerror <- Error.Rate(out_goose$sas[1,] + out_goose$sas[1,], spec = 125, seqlength = 648)


Decile.Plot(out_goose$ownPosCount, seqlength = 648) # plot decile plot
Sliding.Window(out_goose$ownPosCount, seqlength = 648, n = 30) # plot sliding window


### VLF analysis for fishes, N =  2371 barcodes ###

sequences <- fasta.read(file.choose(), 652, pos1 = 1, pos2 = 2)

species.names <- sequences[,2]
specimen.Number <- nrow(sequences)
rownames(sequences) <- species.names
Nuc.count <- count.function(sequences, specimen.Number, seqlength = 652)
frequency.matrix <- ffrequency.matrix.function(Nuc.count, seqlength)
spec.freq <- specimen.frequencies(frequency.matrix, sequences, specimen.Number, species.names, seqlength = 652)
specimen_VLFcount <- VLF.count.spec(spec.freq, p = 0.001, seqlength = 652)
position_VLFcount <- VLF.count.pos(spec.freq, p = 0.001, seqlength = 652)
VLFconvert <- VLF.convert.matrix(sequences, spec.freq, p = 0.001, seqlength = 652)
VLFnuc <- VLF.nucleotides(VLFconvert, sequences, seqlength = 652)
VLFnuc_reduced <- VLF.reduced(VLFnuc, specimen_VLFcount, seqlength = 652)
VLFsingletons <- find.singles(separate(VLFnuc_reduced), seqlength = 652)

sum(specimen_VLFcount) # total number of specimen VLFs - 117
sum(position_VLFcount) # total number of positional VLFs - 117

length(which(specimen_VLFcount != 0)) # total number of specimens with vlfS - 58
length(which(position_VLFcount != 0)) # total number of positions with VLFs - 84

which(specimen_VLFcount != 0) # specimens with VLFs
which(position_VLFcount != 0) # sites containing VLFs

table(specimen_VLFcount) # frequency of specimen VLFs
table(position_VLFcount) # frequency of positional VLFs

single <- Error.Rate(VLFsingletons[1,], VLFsingletons[2,], spec = 2371, seqlength = 652)
shared <- Error.Rate(VLFsingletons[2,], VLFsingletons[1,], spec = 2371, seqlength = 652)
total <- Error.Rate(VLFsingletons[1,] + VLFsingletons[2,], spec = 2371, seqlength = 652)


Decile.Plot(VLFsingletons, seqlength = 652) # plot decile plot
Sliding.Window(VLFsingletons, seqlength = 652) # plot sliding window


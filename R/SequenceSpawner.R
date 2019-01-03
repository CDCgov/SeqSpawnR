#' spawn_sequences
#'
#' @param n The number of sequences to output
#' @param snps The maximum number of snps in each mutant sequence
#' @param seed The original DNA Sequence from which the outputs will be mutated
#'
#' @return A vector of mutant sequences
#' @export
spawn_sequences <- function(n = 20000, snps = 100, seed){
  if(missing(seed)){
    seed <- SeqSpawnR::HXB2
  }
  seednChar <- nchar(seed)

  SampleCodons <- c("GCA", "GCC", "GCG", "GCT", "AAC", "AAT", "GAC", "GAT", "TGC", "TGT", "GAC", "GAT", "GAA", "GAG", "TTC", "TTT", "GGA", "GGC", "GGG", "GGT", "CAC", "CAT", "ATA", "ATC", "ATT" ,"AAA", "AAG" ,"CTA", "CTC", "CTG", "CTT", "TTA", "TTG" , "ATG", "AAC", "AAT", "CCA", "CCC", "CCG", "CCT", "CAA", "CAG", "AGA", "AGG", "CGA", "CGC", "CGG", "CGT", "AGC", "AGT", "TCA", "TCC", "TCG", "TCT", "ACA", "ACC", "ACG", "ACT", "GTA", "GTC", "GTG", "GTT", "TGG", "TAC", "TAT", "CAA", "CAG", "GAA", "GAG" )

  SampleSNPs <- c("A", "C", "G", "T")

  # Set up the vector for the total number of seed variations
  sequences = vector(mode = "character", n)

  #Set up vectors to record the history of the Codon string to search for in the seed sequence, and the replacement
  # condon string
  RandomCodonSetHistory = vector(mode = "character", n)
  ReplacementCodonSetHistory = vector(mode = "character", n)

  #Initialize the seed vector with the initial seed sequence as a seed.
  sequences[1] <- seed

  x <- 1

  while(x < n){
    # The second value can be varied for the total number of possible codons that you want to be considered for consecutive
    # set of condons for variation
    NCodons <- sample(1:10, 1)

    # Of the entire set of codons that represent amino acids randomly select up to NCodons check for existence
    # of the codon string in the current seed under consideration.
    RandomCodonSet <- paste(sample(SampleCodons, NCodons, replace=TRUE), collapse="")
    RandomCodonSetHistory[x] <- RandomCodonSet

    # Do the same and select a codon set of the same length that will replace the searched set.
    ReplacementCodonSet <- paste(sample(SampleCodons, NCodons, replace=TRUE), collapse="")
    ReplacementCodonSetHistory[x] <- ReplacementCodonSet

    # If there is a match between the RandomCodonSet and the seed in processes, set the control flags,
    # and then replace the codon set, update the vector, and increment the count.
    if(RandomCodonSet %>% grep(sequences[x], fixed = TRUE) %>% length > 0){
      Newseed <- str_replace(sequences[sample(1:x, 1)], RandomCodonSet, ReplacementCodonSet)
      sequences[x+1] <- Newseed

      # Add more SNP substitutions randomly across the entire sequence to newly created variant. Replacement numbers are controlled by randomly
      # sampling the AddedSNP number, and randomly picking SNPs to replace.
      for (j in 1:(sample(1:snps, 1))){
        RandomSNP <- sample(SampleSNPs, 1);
        LocOfSNP <- sample(1:seednChar, 1);
        substr(sequences[x+1], LocOfSNP, LocOfSNP) = RandomSNP; # Location of selected SNP to be replaced.
      }
      x <- x+1
    }
  }
  return(sequences)
}

#' Pol Region of the HXB2 Strain of HIV
#'
#' @format A string with 1300 Characters
"HXB2"
#> [1] "HXB2"

#' write_fasta
#'
#' @param sequences A Vector of Sequences
#' @param filename A Filename
#'
#' @export
write_fasta <- function(sequences, filename){
  if(missing(filename)){
    filename <- Sys.time() %>%
      gsub(" ", "-", .) %>%
      gsub(":","-", .) %>%
      paste(".fasta") %>%
      gsub(" ","",.)
  }

  time <- Sys.time()

  n <- length(sequences)

  for(i in 1:n){
    time %>%
      gsub(":", "_", .) %>%
      gsub(" ", "_", .) %>%
      paste(">", .) %>%
      gsub(" ", "", .) %>%
      paste(i) %>%
      gsub(" ", "-", .) %>%
      paste(sequences[i]) %>%
      gsub(" ", "\n", .) %>%
      write(file = filename, append = TRUE)
  }
}

#########################################################################################################
#
#Normalizes PTM to protein level by fitting a global linear model and returning residuals
#MHK June 2020
#
#Fits PTM = beta_0 + beta_1*protein to all matched points in a dataset, and returns residuals as
#protein-corrected PTM values.
#
#The try.all.accession.numbers option will attempt to match all PTM sites with underlying proteins
#using alternative accession numbers.  If there is no initial match in the proteome to a PTM site,
#the function will try all accession numbers associated with the PTM site and pick whichever match
#(if there is one) with the highest protein score.
#
#Requires pacman, cmapR, reshape, ggplot2.  Both proteome and PTM GCTs must have row metadata
#titled 'accession_number' to match features.  If try.all.accession.numbers = TRUE, PTM GCT must have
#accession_numbers' row metadata and proteome GCT must have 'scoreUnique' row metadata.
#
#########################################################################################################

PTM.normalization.lm <- function(
  
  working.directory = "./Input/",
  proteome.gct = 'lscc-v3.2-proteome-ratio-norm-NArm.gct',
  PTM.gct = 'lscc-v3.2-acetylome-ratio-norm-NArm.gct',
  
  try.all.accession.numbers = TRUE,    #will try all PTM accession numbers to find a protein match
  save = FALSE          #save file?
  
)

######################################  No editing necessary below this line  ############################

{
  library('pacman')
  p_load(cmapR)
  p_load(reshape)
  p_load(ggplot2)
  p_load(tidyverse)
  
  
  #set working directory and import GCT files
  setwd(working.directory)
  proteome <- parse_gctx(proteome.gct)
  PTM <- parse_gctx(PTM.gct)

  
  #if PTM accession number does not have match in the proteome, will try to match all accession numbers
  #for that site in the SM output to the proteome. if there are multiple protein matches, it will pick
  #the highest scoring one. replaces value in the 'accession_number' column of the PTM GCT with the new
  #match (duplicates orignal accession number column so those values are saved).
  if(try.all.accession.numbers){
    original_accession_number <- PTM@rdesc$accession_number
    PTM@rdesc <- data.frame(PTM@rdesc, original_accession_number, stringsAsFactors = FALSE)
    PTM@rdesc <- swap.accession.numbers(PTM@rdesc, proteome@rdesc)
  }
  
  
  #fits linear model and returns updated GCT
  PTM.norm <- normalize(PTM, proteome)
  
  
  #writes and returns updated GCT
  if(save){
    file.prefix <- unlist(strsplit(PTM.gct, split = '.gct', fixed = TRUE))[1]
    write_gct(PTM.norm, paste(file.prefix, '_', 'AllAccNum', try.all.accession.numbers,
                              '_proteomenormalized_globalOLS.gct', sep = ''), appenddim = FALSE)
  }
  return(PTM.norm)
  
}



#########################################################################################################
#Matches PTM and proteome by accession number, and tries all possible PTM accession numbers if primary
#accession number doesn't have proteome match

swap.accession.numbers <- function(PTM.rdesc, prot.rdesc){
  prot.rdesc.expanded <- prot.rdesc %>% separate_rows(accession_numbers,sep="\\|")
  
  for(row in c(1:nrow(PTM.rdesc))){
  
    match <- which(prot.rdesc$accession_number == PTM.rdesc$accession_number[row])
    
    if(length(match) == 0){
      #First try to match PTM accession_numbers to protein primary ID
      accession.numbers <- unlist(strsplit(PTM.rdesc$accession_numbers[row],'|', fixed = TRUE))
    
      matches <- prot.rdesc[which(prot.rdesc$accession_number %in% accession.numbers),]
      
      if(nrow(matches) >= 1){
        best <- matches$accession_number[which(matches$scoreUnique == max(matches$scoreUnique, na.rm = TRUE))[1]]
        PTM.rdesc$accession_number[row] <- best
        

      } else {
        #If the above fails, try to match PTM accession_numbers to all Proteome accession_numbers
        matches <- prot.rdesc.expanded[which(prot.rdesc.expanded$accession_numbers %in% accession.numbers),]
        if(nrow(matches) == 1){
          PTM.rdesc$accession_number[row] <- matches$accession_number
        } else if(nrow(matches) > 1){
          best <- matches %>% arrange(desc(scoreUnique)) %>% .$accession_number %>% .[1]
          PTM.rdesc$accession_number[row] <- best
        }
      }
    }
  }
  
  return(PTM.rdesc)
  
}



#########################################################################################################
#Applies linear regression to correct PTM levels for underlying protein levels

normalize <- function(PTM, proteome){
  
  #warn if not all samples are matched
  matched.samples <- intersect(PTM@cid, proteome@cid)
  if(length(matched.samples) != length(PTM@cid)){
    print('WARNING: not all samples in PTM file have matches in proteome.  Unmatched samples removed.')
  }
  
  #create merged data table
  PTM.melt <- melt_gct(PTM)
  prot.melt <- melt_gct(proteome)
  prot.melt.data.only <- data.frame(prot.melt$id.y, prot.melt$accession_number, prot.melt$value)
  colnames(prot.melt.data.only) <- c('id.y', 'accession_number', 'value.prot')
  data <- merge(PTM.melt, prot.melt.data.only, by = c('id.y', 'accession_number'))
  
  #print metrics
  percent <-round(100*nrow(data)/nrow(PTM.melt), digits = 1)
  print(paste(nrow(data), ' points with proteome match out of ', nrow(PTM.melt),
              ' (', percent, '%).', sep = ''))
  
  #fit global model
  print("Fitting model...") 
  model <- lm(value ~ value.prot, data = data)
  residuals <- residuals(model)
  results <- data.frame(data$id.x, data$id.y, residuals)
  colnames(results) <- c('id.x', 'id.y', 'residuals')
  print("Success.")
  print(summary(model))
  
  #compile results into matrix
  results.df <- data.frame(reshape::cast(results, id.x ~ id.y, value.var = residuals))
  results.mat <- data.matrix(results.df[,c(2:ncol(results.df))])
  rownames(results.mat) <- as.character(results.df$id.x)
  
  #reset GCT
  PTM@rdesc <- PTM@rdesc[match(rownames(results.mat), PTM@rdesc$id),]
  PTM@cdesc <- PTM@cdesc[match(colnames(results.mat), PTM@cdesc$id),]
  PTM@rid <- rownames(results.mat)
  PTM@cid <- colnames(results.mat)
  PTM@mat <- results.mat
  
  return(PTM)
}
#This script will calculate the pearson correlation coefficient for all viral regions

#######################################################################################
#######################################################################################

### Sections

# 1. File Set up before running script
# 2. Install &Load required Packages & Import Data
# 3. Run calculation
# 4. Export to csv

#######################################################################################
#######################################################################################

#Section 1: File Set up before running script
#Copy viral regions .fna files to a folder
#Create a single file named X_Eukaryotic_Region.fasta containing a set of large eukaryotic contigs
#Copy X_Eukaryotic_Region.fasta to the viral regions folder 
#The eukaryotic region file has a prefix of X so that it will be the last file in the folder

#######################################################################################
#######################################################################################

# 2. Install &Load required Packages & Import Data
#    The TNF_Set is essentially a matrix containing the TNF for all viral regions

install.packages("biostrings")

library(Biostrings)
 vRegions<-list.files(path="C:\\Users\\USERNAME\\FOLDERNAME", full.names=TRUE)
 files <- open_input_files(vRegions)
 DNA_Set<-readDNAStringSet(files)
 TNF_Set <- oligonucleotideFrequency(DNA_Set,4,step=4)
 TNF_Set[1,]                                                
 print(length(files))
 size <- length(files)
 

# cor is a built in function for finding pearson correlation
# It could be used in the main loop instead of the long pearson calculation
# Its included here just as a sanity check to confirm the long pearson calc is correct
 res <- cor(TNF_Set[1,], TNF_Set[size, ], method = "pearson")
 print(res)                                           

#######################################################################################
#######################################################################################
# 3. Run calculation
#     
pearson<-0
for(i in 1:(size-1)){
        pearson[i] <-   sum (   (TNF_Set[i,] - mean(TNF_Set[i,])) * (TNF_Set[size,] - mean(TNF_Set[size,]))  ) / 
        sqrt (   sum ( (TNF_Set[i,] - mean(TNF_Set[i,])) * (TNF_Set[i,] - mean(TNF_Set[i,])) ) *  
                         sum ( (TNF_Set[size,] - mean(TNF_Set[size,])) * (TNF_Set[size,] - mean(TNF_Set[size,])) )    ) 
       
        con <- file(vRegions[i],"r")
        first_line <- readLines(con,n=1)
        close(con)

        names(pearson)[i]<-paste0(first_line)
       # print(pearson[i]) 
        
}

#######################################################################################
#######################################################################################
# 4. Export sorted results to csv

sPearson <-sort(pearson)
print(pearson)
write.table(sPearson, file = "C:\\Users\\USERNAME\\FOLDERNAME\\Output\\output.csv", col.names=FALSE)

#######################################################################################
#######################################################################################







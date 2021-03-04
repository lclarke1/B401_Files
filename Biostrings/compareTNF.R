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

#######################################################################################
#######################################################################################

# 2. Install &Load required Packages & Import Data
#    The TetSocTNF is essentially a matrix containing the TNF for all viral regions

install.packages("biostrings")

library(Biostrings)
 vRegions<-list.files(path="C:\\Users\\USERNAME\\FOLDERNAME", full.names=TRUE)
 files <- open_input_files(vRegions)
 TetSoc<-readDNAStringSet(files)
 TetSocTNF <- oligonucleotideFrequency(TetSoc,4, step=300)
 TetSocTNF[1,]                                                

#######################################################################################
#######################################################################################
# 3. Run calculation
#    Note: You will need to adjust the calculation range.
#          The upper range of the for loop will be the number of viral regions in the folder
#          The index here for file 93 is referring to the eukaryotic regions file
#          Example
#          If you have 10 viral regions, and the eukarotic region file is the 11th file
#          The for loop range will be 1:10 and all entries that reference region 93 become 11. 

pearson<-0
for(i in 1:92){
        pearson[i] <-   sum (   (TetSocTNF[i,] - mean(TetSocTNF[i,])) * (TetSocTNF[93,] - mean(TetSocTNF[93,]))  ) / 
        sqrt (   sum ( (TetSocTNF[i,] - mean(TetSocTNF[i,])) * (TetSocTNF[i,] - mean(TetSocTNF[i,])) ) *  
                         sum ( (TetSocTNF[93,] - mean(TetSocTNF[93,])) * (TetSocTNF[93,] - mean(TetSocTNF[93,])) )    ) 
       
        con <- file(vRegions[i],"r")
        first_line <- readLines(con,n=1)
        close(con)

        names(pearson)[i]<-paste0(first_line)
        print(pearson[i]) 
        
}

#######################################################################################
#######################################################################################
# 4. Export sorted results to csv

sPearson <-sort(pearson)
print(sPearson)
write.table(sPearson, file = "C:\\Users\\USERNAME\\FOLDERNAME\\Output\\output.csv", col.names=FALSE)

#######################################################################################
#######################################################################################







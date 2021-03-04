 library(Biostrings)
 vRegions<-list.files(path="C:\\Users\\USERNAME\\FOLDERNAME", full.names=TRUE)
 files <- open_input_files(vRegions)
 TetSoc<-readDNAStringSet(files)
 TetSocTNF <- oligonucleotideFrequency(TetSoc,4, step=300)
 TetSocTNF[1,]                                                

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

sPearson <-sort(pearson)
print(sPearson)
write.table(sPearson, file = "C:\\Users\\USERNAME\\FOLDERNAME\\Output\\output.csv", col.names=FALSE)



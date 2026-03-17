#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# R_utils.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Utils scripts for R
#
# Author: Jurriaan Janssen (j.janssen.1@erasmusmc.nl)
#
# condaenv: 
# Usage: 
#
# TODO:
# 1) 
#
# History:
#  17-03-2026: File creation
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1.0  Specify functions
#-------------------------------------------------------------------------------
get_color_palette <- function(pal='standard'){
    if(pal=='standard'){
        # Tableau 10
        color_pal <-  c("#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F","#EDC948","#B07AA1","#FF9DA7","#9C755F","#BAB0AC")
        names(color_pal) <- c("Blue",  "Orange",  "Red",  "Light Teal","Green","Yellow", "Purple", "Pink",   "Brown", "Light Gray")
    }else if(pal == 'paired'){
        # tableau 20
        color_pal <- c("#4E79A7","#A0CBE8","#F28E2B","#FFBE7D","#59A14F","#8CD17D","#B6992D","#F1CE63","#499894","#86BCB6","#E15759","#FF9D9A","#79706E","#BAB0AC","#D37295","#FABFD2","#B07AA1","#D4A6C8","#9D7660","#D7B5A6")
        names(color_pal) <- c("Blue","Light Blue" ,"Orange","Light Orange","Green","Light Green","Yellow-Green","Yellow","Teal","Light Teal", "Red","Pink","Dark Gray","Light Gray","Pink","Light Pink","Purple","Light Purple","Brown","Light Orange")
    }
    return(color_pal)
}


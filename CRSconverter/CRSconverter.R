#!/bin/Rscript
# GRID CREATION
args = commandArgs(trailingOnly=TRUE)

library(sf)
library(tidyr)


save_as_pdf_fun <- function(){

    pdf(file = "output.pdf")
    plot(st_geometry(Transformed_spacial_coordinates))
    dev.off()
}


save_as_image_fun <- function(format){

    call(format, filename = paste("output.",format), sep = "", collapse = NULL)
    print(call(format, filename = paste("output.",format)), sep = "", collapse = NULL)
    plot(st_geometry(Transformed_spacial_coordinates))
    dev.off()
}


save_as_shp_fun <- function(){
    write_sf(Transformed_spacial_coordinates, "output.shp")
}




if (length(args)<1){stop('please provide spacial coordinates files(.shp files)')
}else{
    Spacial_coordinates_files <- read_sf(as.character(args[1]), layer = 'shapefile')
    projection <-    paste('+proj='
                          , as.character(args[2])
                          , ' +lat_0='
                          , as.character(args[3])
                          , ' +lon_0='
                          , as.character(args[4])
                          , ' +x_0='
                          , as.character(args[5])
                          , ' +y_0='
                          , as.character(args[6])
                          , ' +ellps='
                          , as.character(args[7])
                          , ' +datum='
                          , as.character(args[7])
                          , ' +units='
                          , as.character(args[8])
                          , ' +no_defs', sep = "", collapse = NULL) #presonalisation du systeme de références

    Transformed_spacial_coordinates <- st_transform(Spacial_coordinates_files, crs = projection)
    
    for (a in 9:length(args)) {
        
        if (as.character(args[a]) == "pdf") {
        save_as_pdf_fun()
        next

        }else if (as.character(args[a]) == "shp") {
        save_as_shp_fun()
        next

        }else { save_as_image_fun(as.character(args[a])) }
    }}






#write.table(Transformed_spacial_coordinates, file = "output.tabular", sep="\t", row.names=FALSE)
    #enregistrement dans quoi ?
    
#this tool convert shp files to the right format for whatever you like ~a phylodiversity analysis

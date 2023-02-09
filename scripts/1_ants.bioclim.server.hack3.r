#################
require(AntWeb)
require(raster)
require(assertthat)
require(rjson)
require(plyr)

#allows increased limit - function originally from AntWeb
aw_data_hack<-function (genus = NULL, species = NULL, scientific_name = NULL, 
    georeferenced = NULL, min_elevation = NULL, max_elevation = NULL, 
    type = NULL, habitat = NULL, country = NULL, min_date = NULL, 
    max_date = NULL, bbox = NULL, limit = NULL, offset = NULL, 
    quiet = FALSE) 
{
    require(AntWeb)
    require(assertthat)
    require(rjson)
    main_args <- AntWeb:::z_compact(as.list(c(scientific_name, genus, 
        type, habitat, bbox)))
    date_args <- AntWeb:::z_compact(as.list(c(min_date, max_date)))
    elev_args <- AntWeb:::z_compact(as.list(c(min_elevation, max_elevation)))
    arg_lengths <- c(length(main_args), length(date_args), length(elev_args))
    assert_that(any(arg_lengths) > 0)
    decimal_latitude <- NA
    decimal_longitude <- NA
    if (!is.null(scientific_name)) {
        genus <- strsplit(scientific_name, " ")[[1]][1]
        species <- strsplit(scientific_name, " ")[[1]][2]
    }
    base_url <- "http://www.antweb.org/api/v2/"
    original_limit <- limit
    args <- AntWeb:::z_compact(as.list(c(genus = genus, species = species, 
        bbox = bbox, min_elevation = min_elevation, max_elevation = max_elevation, 
        habitat = habitat, country = country, type = type, min_date = min_date, 
        max_date = max_date, limit = 1, offset = offset, georeferenced = georeferenced)))
    results <- httr:::GET(base_url, query = args)
    httr:::warn_for_status(results)
    data <- rjson:::fromJSON(httr:::content(results, "text"))
    data <- AntWeb:::z_compact(data)
    if (data$count > 100000 & is.null(limit)) {
        args$limit <- 100000
        results <- httr:::GET(base_url, query = args)
        if (!quiet) 
            message(sprintf("Query contains %s results. First 1000 retrieved. Use the offset argument to retrieve more \n", 
                data$count))
    }
    else {
        args$limit <- original_limit
        results <- httr:::GET(base_url, query = args)
    }
    data <- rjson:::fromJSON(httr:::content(results, "text"))
    data <- AntWeb:::z_compact(data)
    if (identical(data$specimens$empty_set, "No records found.")) {
        NULL
    }
    else {
        if (!quiet) 
            message(sprintf("%s results available for query.", 
                data$count))
        data_df <- lapply(data$specimens, function(x) {
            x$images <- NULL
            df <- data.frame(t(unlist(x)), stringsAsFactors = FALSE)
            df
        })
        final_df <- data.frame(do.call(plyr:::rbind.fill, data_df))
        names(final_df)[grep("geojson.coord1", names(final_df))] <- "decimal_latitude"
        names(final_df)[grep("geojson.coord2", names(final_df))] <- "decimal_longitude"
        final_df$decimalLatitude <- NULL
        final_df$decimalLongitude <- NULL
        final_df$minimumElevationInMeters <- as.numeric(final_df$minimumElevationInMeters)
        final_results <- list(count = data$count, call = args, 
            data = final_df)
        class(final_results) <- "antweb"
        final_results
    }
}

antweb.georef.grabber<-function(master.list,group.size){
	master.list$AntWebPresent<-NA
	master.list$georeferenced<-NA
	n.reps<-ceiling(nrow(master.list)/group.size)
	ls<-getData(name="worldclim",var="bio", res=2.5, download=TRUE)
	for(n in 1:n.reps){
		ant.dat<-NULL
		attempt<-1
		no.samp<-0
		while(no.samp<(n*group.size) && attempt<=5){
			if(attempt>1){
				Sys.sleep(runif(1,min=1,max=10))
			}
			for(x in (1+(group.size*(n-1))):(n*group.size)){
				print(paste("Searching for",master.list$Taxon[x],": Taxon", x, sep=" "))
				ad<-NULL
				try(ad<-aw_data_hack(scientific_name=gsub("_"," ",master.list$Taxon[x])))
				if(!is.null(ad)){
					master.list$AntWebPresent[x]<-"Present"
					print(paste("Searching for georeferenced",master.list$Taxon[x],": Taxon", x, sep=" "))
					ant.dat<-NULL
					try(ant.dat<-aw_data_hack(scientific_name=gsub("_"," ",master.list$Taxon[x]),georeferenced=TRUE))
					if(!is.null(ant.dat)){
						master.list$georeferenced[x]<-"Present"
						ant.dat$data$decimal_longitude<-as.numeric(ant.dat$data$decimal_longitude)
						ant.dat$data$decimal_latitude<-as.numeric(ant.dat$data$decimal_latitude)
						ant.dat$data[,paste("bio",1:19,sep="")]<-NA
						ant.dat$data[,paste("bio",1:19,sep="")]<-extract(ls,ant.dat$data[,c("decimal_longitude","decimal_latitude")])
						write.csv(ant.dat$data,file=paste("/home/mpnelsen/antweb_bioclim_vars_06sep2017/individual_files/",paste(master.list$Taxon[x],"bioclim_vars.csv",sep="_"),sep=""),row.names=FALSE)
					}	
				}	
				no.samp<-x+1
			}
			attempt<-attempt+1
			#print(attempt)
			#print(no.samp)
			#print(no.samp<(n*group.size) && attempt<=5)
		}
	}
return(master.list)
}

sink(file="/home/mpnelsen/antweb_bioclim_vars_06sep2017/screen_output.txt",append=TRUE,type="output")
ants<-read.csv(file="/home/mpnelsen/antweb_bioclim_vars_06sep2017/AntWiki_Valid_Species_w_fossils_mymods_EXTANT_ONLY_WEXTRAGENUSPLACEH.csv",stringsAsFactors=FALSE)
new.ants<-antweb.georef.grabber(ants,100)
head(new.ants)
write.csv(new.ants,file="/home/mpnelsen/antweb_bioclim_vars_06sep2017/AntWiki_Valid_Species_w_fossils_mymods_EXTANT_ONLY_WEXTRAGENUSPLACEH_w_awinfo.csv",row.names=FALSE)
sink()



#Introduced Ants copied from: https://www.antweb.org/taxonomicPage.do?rank=species&project=introducedants on 28 December 2016 from AntWeb v6.13.7 cleaned and saved as antweb.introduced.28dec2016.clean.csv
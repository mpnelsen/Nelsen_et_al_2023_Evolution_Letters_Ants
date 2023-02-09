#modified from rangeBuilder
standardizeCountryHack<-function (country, fuzzyDist = 1, nthreads = 1) 
{
	countryList<-rangeBuilder:::countryList
    if (any(is.na(country))) {
        country[which(is.na(country))] <- ""
    }
    country <- toupper(country)
    country <- gsub("_|-", " ", country)
    country <- gsub("\\\\.", "", country)
    country <- stringi::stri_trans_general(country, "Latin-ASCII")
    country <- gsub("(^|\\\\s)ST\\\\.?\\\\s", "SAINT ", country)
    country <- gsub("\\\\?|\\\\[|\\\\]", "", country)
    country <- gsub("\\\\/", "", country)
    country <- gsub("\\\\s+", " ", country)
    country <- trim(country)
    res <- character(length(country))
    matchCountry <- function(val, countryList, fuzzyDist = fuzzyDist) {
        if (nchar(val) == 2) {
            if (val %in% isoLookup[, 2]) {
                ind <- which(names(countryList) == isoLookup[which(isoLookup[, 
                  2] == val), 1])
            }
            else {
                ind <- NULL
            }
        }
        if (nchar(val) == 3) {
            if (val %in% isoLookup[, 3]) {
                ind <- which(names(countryList) == isoLookup[which(isoLookup[, 
                  3] == val), 1])
            }
            else {
                ind <- NULL
            }
        }
        if (nchar(val) > 3) {
            ind <- which(sapply(countryList, function(y) val %in% 
                y) == TRUE)
            if (length(ind) == 0) {
                d <- adist(val, names(countryList))
                mind <- which.min(d)
                if (d[mind] <= fuzzyDist) {
                  ind <- mind
                }
                else {
                  d <- sapply(countryList, function(y) adist(val, 
                    y))
                  d <- sapply(d, min)
                  mind <- which.min(d)
                  if (d[mind] <= fuzzyDist) {
                    ind <- mind
                  }
                }
            }
        }
        if (nchar(val) < 2) {
            ind <- NULL
        }
        if (length(ind) == 0) {
            return("")
        }
        else {
            return(names(countryList)[ind])
        }
    }
    uniqueCountry <- unique(country)
    if (nthreads > 1) {
        cl <- parallel::makePSOCKcluster(nthreads)
        parallel::clusterExport(cl = cl, varlist = c("uniqueCountry", 
            "countryList"), envir = environment())
        uniqueRes <- sapply(uniqueCountry, function(x) {
            return(matchCountry(x, countryList, fuzzyDist = fuzzyDist))
        }, simplify = TRUE, USE.NAMES = FALSE, cl = cl)
        parallel::stopCluster(cl)
    }
    else {
        uniqueRes <- sapply(uniqueCountry, function(x) {
            return(matchCountry(x, countryList, fuzzyDist = fuzzyDist))
        }, simplify = TRUE, USE.NAMES = FALSE)
    }
    for (i in 1:length(uniqueCountry)) {
        res[which(country == uniqueCountry[i])] <- uniqueRes[i]
    }
    rm(countryList)
    rm(uniqueCountry)
    rm(uniqueRes)
    rm(country)
    return(res)
}



vars<-c("biomes","anntotprecip","avgannrh","elevation","evapotrans","npp","pevapotrans")

file.path<-"~/Desktop/ants_ecogeography/antweb_bioclim_vars_06sep2017/individual_files_w_addl_data/"
fls<-list.files(path=file.path)
fls<-gsub("_bioclim_vars.csv","",fls)


#get list of countries
countries<-NA
for(x in 1:length(fls)){
	df<-read.csv(file=paste(file.path,fls[x],"_bioclim_vars.csv",sep=""),stringsAsFactors=FALSE)
	if(nrow(df)>=99999){
		print(paste("Limit potentially reached for ",fls[x],sep=""))
	}
	tax.in.file<-paste(df$genus[1],df$specificEpithet[1],sep="_")
	if(tax.in.file!=fls[x]){
		print(paste("Problem with ",fls[x]," file - contains info for ",tax.in.file,sep=""))
	}	
	if(tax.in.file==fls[x]){
		countries<-c(countries,unique(df$country))
	}
countries<-(sort(unique(countries)))
}	

library(rangeBuilder)
#Used rangeBuilder 1.4
standardized.countries<-standardizeCountry(countries)


file.path.mod<-"~/Desktop/ants_ecogeography/antweb_bioclim_vars_06sep2017/individual_files_w_addl_data_countrymods/"
problems<-data.frame(matrix(nrow=0,ncol=5))
colnames(problems)<-c("Taxon","decimal_longitude","decimal_latitude","country_stand","closest_country")
#zz<-0

#require(geonames)
#require(countrycode)
#options(geonamesUsername="mnelsen")
require(maps)
require(rangeBuilder)

#This updates countries
for(x in 1:length(fls)){
	df<-read.csv(file=paste(file.path,fls[x],"_bioclim_vars.csv",sep=""),stringsAsFactors=FALSE)
	df[,c("country_mod","country_stand","closest_country","decimal_longitude_corr","decimal_latitude_corr","coords_flipped")]<-NA
	df[,"coords_flipped"]<-"No"
	df$Land<-NA
	df$Keep<-"No"
	if(nrow(df)>=99999){
		print(paste("Limit potentially reached for ",fls[x],sep=""))
	}
	tax.in.file<-paste(df$genus[1],df$specificEpithet[1],sep="_")
	if(tax.in.file!=fls[x]){
		print(paste("Problem with ",fls[x]," file - contains info for ",tax.in.file,sep=""))
	}	
	if(tax.in.file==fls[x]){
		df$country_mod<-df$country
		df$country_mod<-gsub("C te d Ivoire","Cote d'Ivoire",df$country_mod)
		df$country_mod<-gsub("Democratic Republic of Congo","Democratic Republic of the Congo",df$country_mod)
		#df$country_mod[df$country_mod=="Congo"]<-"Democratic Republic of the Congo"
		df$country_mod<-gsub("ECUADOR","Ecuador",df$country_mod)
		df$country_mod<-gsub("Federated States of Micronesia","Micronesia",df$country_mod)
		df$country_mod<-gsub("Polinesia","Polynesia",df$country_mod)
		df$country_mod<-gsub("Rhodesia","Zimbabwe",df$country_mod)
		#df$country_mod<-gsub("Samoa","American Samoa",df$country_mod)
		df$country_mod<-gsub("Siam Thailand","Thailand",df$country_mod)
		df$country_mod<-gsub("United States of America","United States",df$country_mod)
		df$country_mod<-gsub("Viet Nam","Vietnam",df$country_mod)
		#df$country_mod<-gsub("DEMOCRATIC REPUBLIC OF THE CONGO","CONGO",df$country_mod)
		df$country_stand<-standardizeCountryHack(df$country_mod)
		df$country_stand[df$country_stand==""]<-NA
		for(p in 1:nrow(df)){
			#print(paste(df$decimal_latitude[p],df$decimal_longitude[p],sep=" "))
			#clos.full<-GNcountryCode(df$decimal_latitude[p],df$decimal_longitude[p])
			#clos<-toupper(countrycode(clos.full$countryCode,"iso2c","country.name"))
			#clos<-closestCountry(df[p,c("decimal_longitude","decimal_latitude")])
			clos<-strsplit(toupper(map.where(database="world",x=df$decimal_longitude[p],y=df$decimal_latitude[p])),":")[[1]][1]
			clos[clos=="USA"]<-"UNITED STATES"
			clos[clos=="TRINIDAD"]<-"TRINIDAD AND TOBAGO"
			clos[clos=="TOBAGO"]<-"TRINIDAD AND TOBAGO"
			#clos[clos=="SAMOA"]<-"AMERICAN SAMOA"
			clos[clos=="SAINT VINCENT"]<-"SAINT VINCENT AND THE GRENADINES"
			clos[clos=="GRENADINES"]<-"SAINT VINCENT AND THE GRENADINES"
			clos[clos=="IVORY COAST"]<-"COTE D'IVOIRE"
			#clos[clos=="REPUBLIC OF CONGO"]<-"DEMOCRATIC REPUBLIC OF THE CONGO"
			clos[clos=="UK"]<-"UNITED KINGDOM"
	
			#Exceptions allowed if:
			clos[clos=="CANARY ISLANDS" & df[p,"country_stand"]=="SPAIN"]<-"SPAIN"
			clos[clos=="MEXICO" & df[p,"country_stand"]=="UNITED STATES"]<-"UNITED STATES"
			clos[clos=="HONDURAS" & df[p,"country_stand"]=="NICARAGUA"]<-"NICARAGUA"
			clos[clos=="INDONESIA" & df[p,"country_stand"]=="MALAYSIA"]<-"MALAYSIA"
			clos[clos=="GUATEMALA" & df[p,"country_stand"]=="MEXICO"]<-"MEXICO"
			clos[clos=="NORTHERN MARIANA ISLANDS" & df[p,"country_stand"]=="MICRONESIA"]<-"MICRONESIA"
			clos[clos=="NIUE" & df[p,"country_stand"]=="COOK ISLANDS"]<-"COOK ISLANDS"
			clos[clos=="CHRISTMAS ISLAND" & df[p,"country_stand"]=="AUSTRALIA"]<-"AUSTRALIA"
			clos[clos=="COSTA RICA" & df[p,"country_stand"]=="NICARAGUA"]<-"NICARAGUA"
			clos[clos=="UNITED STATES" & df[p,"country_stand"]=="CANADA"]<-"CANADA"
			clos[clos=="UNITED ARAB EMIRATES" & df[p,"country_stand"]=="OMAN"]<-"OMAN"
			clos[clos=="FRANCE" & df[p,"country_stand"]=="BELGIUM"]<-"BELGIUM"
			clos[clos=="PUERTO RICO" & df[p,"country_stand"]=="UNITED STATES"]<-"UNITED STATES"
			clos[clos=="TANZANIA" & df[p,"country_stand"]=="ZAMBIA"]<-"ZAMBIA"
			clos[clos=="MOROCCO" & df[p,"country_stand"]=="WESTERN SAHARA"]<-"WESTERN SAHARA"
			clos[clos=="CANADA" & df[p,"country_stand"]=="UNITED STATES"]<-"UNITED STATES"
			clos[clos=="CHINA" & df[p,"country_stand"]=="MONGOLIA"]<-"MONGOLIA"
			clos[clos=="LUXEMBOURG" & df[p,"country_stand"]=="BELGIUM"]<-"BELGIUM"
			clos[clos=="GERMANY" & df[p,"country_stand"]=="BELGIUM"]<-"BELGIUM"
			clos[clos=="SWITZERLAND" & df[p,"country_stand"]=="ITALY"]<-"ITALY"
			clos[clos=="BRAZIL" & df[p,"country_stand"]=="PERU"]<-"PERU"
			clos[clos=="CROATIA" & df[p,"country_stand"]=="SLOVENIA"]<-"SLOVENIA"
			clos[clos=="ANGOLA" & df[p,"country_stand"]=="DEMOCRATIC REPUBLIC OF THE CONGO"]<-"DEMOCRATIC REPUBLIC OF THE CONGO"
			clos[clos=="OMAN" & df[p,"country_stand"]=="UNITED ARAB EMIRATES"]<-"UNITED ARAB EMIRATES"
			clos[clos=="PANAMA" & df[p,"country_stand"]=="COSTA RICA"]<-"COSTA RICA"
			clos[clos=="ARGENTINA" & df[p,"country_stand"]=="BRAZIL"]<-"BRAZIL"
			clos[clos=="PARAGUAY" & df[p,"country_stand"]=="ARGENTINA"]<-"ARGENTINA"
			clos[clos=="UNITED STATES" & df[p,"country_stand"]=="MEXICO"]<-"MEXICO"
			clos[clos=="SURINAME" & df[p,"country_stand"]=="GUYANA"]<-"GUYANA"
			clos[clos=="VENEZUELA" & df[p,"country_stand"]=="COLOMBIA"]<-"COLOMBIA"
			clos[clos=="CAYMAN ISLANDS" & df[p,"country_stand"]=="UNITED KINGDOM"]<-"UNITED KINGDOM"
			clos[clos=="PALAU" & df[p,"country_stand"]=="MICRONESIA"]<-"MICRONESIA"
			clos[clos=="COLOMBIA" & df[p,"country_stand"]=="PERU"]<-"PERU"
			clos[clos=="BRAZIL" & df[p,"country_stand"]=="BOLIVIA"]<-"BOLIVIA"
			clos[clos=="HAITI" & df[p,"country_stand"]=="UNITED STATES"]<-"UNITED STATES"
			clos[clos=="UNITED KINGDOM" & df[p,"country_stand"]=="SCOTLAND"]<-"SCOTLAND"
			clos[clos=="NETHERLANDS" & df[p,"country_stand"]=="BELGIUM"]<-"BELGIUM"
			clos[clos=="UNITED KINGDOM" & df[p,"country_stand"]=="ENGLAND"]<-"ENGLAND"
			clos[clos=="GERMANY" & df[p,"country_stand"]=="CZECH REPUBLIC"]<-"CZECH REPUBLIC"
			clos[clos=="TURKMENISTAN" & df[p,"country_stand"]=="IRAN"]<-"IRAN"
			clos[clos=="SOUTH AFRICA" & df[p,"country_stand"]=="BOTSWANA"]<-"BOTSWANA"
			clos[clos=="ETHIOPIA" & df[p,"country_stand"]=="SOMALIA"]<-"SOMALIA"
			clos[clos=="MADEIRA ISLANDS" & df[p,"country_stand"]=="PORTUGAL"]<-"PORTUGAL"
			clos[clos=="URUGUAY" & df[p,"country_stand"]=="ARGENTINA"]<-"ARGENTINA"
			clos[clos=="ARGENTINA" & df[p,"country_stand"]=="PARAGUAY"]<-"PARAGUAY"
			clos[clos=="AZORES" & df[p,"country_stand"]=="PORTUGAL"]<-"PORTUGAL"
			clos[clos=="VIRGIN ISLANDS, US" & df[p,"country_stand"]=="UNITED STATES VIRGIN ISLANDS"]<-"UNITED STATES VIRGIN ISLANDS"
			clos[clos=="GUAM" & df[p,"country_stand"]=="MICRONESIA"]<-"MICRONESIA"
			clos[clos=="CAMBODIA" & df[p,"country_stand"]=="VIETNAM"]<-"VIETNAM"
			clos[clos=="NORFOLK ISLAND" & df[p,"country_stand"]=="AUSTRALIA"]<-"AUSTRALIA"
			clos[clos=="SOUTH AFRICA" & df[p,"country_stand"]=="ZIMBABWE"]<-"ZIMBABWE"
			clos[clos=="SOUTH SUDAN" & df[p,"country_stand"]=="SUDAN"]<-"SUDAN"
			clos[clos=="BOLIVIA" & df[p,"country_stand"]=="PERU"]<-"PERU"
			clos[clos=="MALAYSIA" & df[p,"country_stand"]=="INDONESIA"]<-"INDONESIA"
			clos[clos=="INDONESIA" & df[p,"country_stand"]=="PAPUA NEW GUINEA"]<-"PAPUA NEW GUINEA"
			clos[clos=="EL SALVADOR" & df[p,"country_stand"]=="GUATEMALA"]<-"GUATEMALA"
			clos[clos=="BRUNEI" & df[p,"country_stand"]=="MALAYSIA"]<-"MALAYSIA"
			clos[clos=="MONTENEGRO" & df[p,"country_stand"]=="CROATIA"]<-"CROATIA"
			clos[clos=="BOLIVIA" & df[p,"country_stand"]=="ARGENTINA"]<-"ARGENTINA"
			clos[clos=="MYANMAR" & df[p,"country_stand"]=="THAILAND"]<-"THAILAND"
			clos[clos=="MEXICO" & df[p,"country_stand"]=="GUATEMALA"]<-"GUATEMALA"
			clos[clos=="PERU" & df[p,"country_stand"]=="ECUADOR"]<-"ECUADOR"
			clos[clos=="COLOMBIA" & df[p,"country_stand"]=="BRAZIL"]<-"BRAZIL"
			clos[clos=="BRAZIL" & df[p,"country_stand"]=="COLOMBIA"]<-"COLOMBIA"
			clos[clos=="LESOTHO" & df[p,"country_stand"]=="SOUTH AFRICA"]<-"SOUTH AFRICA"
			clos[clos=="ZAMBIA" & df[p,"country_stand"]=="ZIMBABWE"]<-"ZIMBABWE"
			clos[clos=="VIRGIN ISLANDS, BRITISH" & df[p,"country_stand"]=="BRITISH VIRGIN ISLANDS"]<-"BRITISH VIRGIN ISLANDS"
			clos[clos=="SAMOA" & df[p,"country_stand"]=="AMERICAN SAMOA"]<-"AMERICAN SAMOA"
			clos[clos=="REPUBLIC OF CONGO" & df[p,"country_stand"]=="CONGO"]<-"CONGO"
			clos[clos=="DEMOCRATIC REPUBLIC OF THE CONGO" & df[p,"country_stand"]=="CONGO"]<-"CONGO"
			if(df[p,"country_stand"] %in% clos){
				df[p,c("decimal_longitude_corr","decimal_latitude_corr")]<-df[p,c("decimal_longitude","decimal_latitude")]
				df[p,"closest_country"]<-df[p,"country_stand"]
				df$Land[p]<-filterByLand(df[p,c("decimal_longitude_corr","decimal_latitude_corr")])
				if(isTRUE(df$Land[p]==TRUE & df$country_stand[p]==df$closest_country[p])){
					df$Keep[p]<-"Yes"
				}
				#if(df$Keep[p]=="No"){
					#print(paste(df$Taxon[p],df$decimal_longitude[p],df$decimal_latitude[p],"listed in and matches", df$country_stand[p],"but not over land",sep=" "))
				#}
			}
			if(!df[p,"country_stand"] %in% clos){
				flip.df<-as.data.frame(matrix(data=NA,nrow=4,ncol=3))
				colnames(flip.df)<-c("decimal_latitude","decimal_longitude","Co")
				flip.df$decimal_latitude<-df[p,c("decimal_latitude")]*c(1,1,-1,-1)
				flip.df$decimal_longitude<-df[p,c("decimal_longitude")]*c(1,-1,1,-1)
				flip.df$Co<-toupper(map.where(database="world",x=flip.df$decimal_longitude,y=flip.df$decimal_latitude))
				flip.df$Co<-sapply(strsplit(flip.df$Co,":"),"[",1)
				flip.df$Co[flip.df$Co=="USA"]<-"UNITED STATES"
				flip.df$Co[flip.df$Co=="TRINIDAD"]<-"TRINIDAD AND TOBAGO"
				flip.df$Co[flip.df$Co=="TOBAGO"]<-"TRINIDAD AND TOBAGO"
				#flip.df$Co[flip.df$Co=="SAMOA"]<-"AMERICAN SAMOA"
				#flip.df$Co[flip.df$Co=="AMERICAN AMERICAN SAMOA"]<-"AMERICAN SAMOA"
				flip.df$Co[flip.df$Co=="SAINT VINCENT"]<-"SAINT VINCENT AND THE GRENADINES"
				flip.df$Co[flip.df$Co=="GRENADINES"]<-"SAINT VINCENT AND THE GRENADINES"
				flip.df$Co[flip.df$Co=="IVORY COAST"]<-"COTE D'IVOIRE"
				#flip.df$Co[flip.df$Co=="REPUBLIC OF CONGO"]<-"DEMOCRATIC REPUBLIC OF THE CONGO"
				#flip.df$Co[flip.df$Co=="CANARY ISLANDS"]<-"SPAIN"
				flip.df$Co[flip.df$Co=="UK"]<-"UNITED KINGDOM"

				#Exceptions allowed if:
				flip.df$Co[flip.df$Co=="CANARY ISLANDS" & df[p,"country_stand"]=="SPAIN"]<-"SPAIN"
				flip.df$Co[flip.df$Co=="MEXICO" & df[p,"country_stand"]=="UNITED STATES"]<-"UNITED STATES"
				flip.df$Co[flip.df$Co=="HONDURAS" & df[p,"country_stand"]=="NICARAGUA"]<-"NICARAGUA"
				flip.df$Co[flip.df$Co=="INDONESIA" & df[p,"country_stand"]=="MALAYSIA"]<-"MALAYSIA"
				flip.df$Co[flip.df$Co=="GUATEMALA" & df[p,"country_stand"]=="MEXICO"]<-"MEXICO"
				flip.df$Co[flip.df$Co=="NORTHERN MARIANA ISLANDS" & df[p,"country_stand"]=="MICRONESIA"]<-"MICRONESIA"
				flip.df$Co[flip.df$Co=="NIUE" & df[p,"country_stand"]=="COOK ISLANDS"]<-"COOK ISLANDS"
				flip.df$Co[flip.df$Co=="CHRISTMAS ISLAND" & df[p,"country_stand"]=="AUSTRALIA"]<-"AUSTRALIA"
				flip.df$Co[flip.df$Co=="COSTA RICA" & df[p,"country_stand"]=="NICARAGUA"]<-"NICARAGUA"
				flip.df$Co[flip.df$Co=="UNITED STATES" & df[p,"country_stand"]=="CANADA"]<-"CANADA"
				flip.df$Co[flip.df$Co=="UNITED ARAB EMIRATES" & df[p,"country_stand"]=="OMAN"]<-"OMAN"
				flip.df$Co[flip.df$Co=="FRANCE" & df[p,"country_stand"]=="BELGIUM"]<-"BELGIUM"
				flip.df$Co[flip.df$Co=="PUERTO RICO" & df[p,"country_stand"]=="UNITED STATES"]<-"UNITED STATES"
				flip.df$Co[flip.df$Co=="TANZANIA" & df[p,"country_stand"]=="ZAMBIA"]<-"ZAMBIA"
				flip.df$Co[flip.df$Co=="MOROCCO" & df[p,"country_stand"]=="WESTERN SAHARA"]<-"WESTERN SAHARA"
				flip.df$Co[flip.df$Co=="CANADA" & df[p,"country_stand"]=="UNITED STATES"]<-"UNITED STATES"
				flip.df$Co[flip.df$Co=="CHINA" & df[p,"country_stand"]=="MONGOLIA"]<-"MONGOLIA"
				flip.df$Co[flip.df$Co=="LUXEMBOURG" & df[p,"country_stand"]=="BELGIUM"]<-"BELGIUM"
				flip.df$Co[flip.df$Co=="GERMANY" & df[p,"country_stand"]=="BELGIUM"]<-"BELGIUM"
				flip.df$Co[flip.df$Co=="SWITZERLAND" & df[p,"country_stand"]=="ITALY"]<-"ITALY"
				flip.df$Co[flip.df$Co=="BRAZIL" & df[p,"country_stand"]=="PERU"]<-"PERU"
				flip.df$Co[flip.df$Co=="CROATIA" & df[p,"country_stand"]=="SLOVENIA"]<-"SLOVENIA"
				flip.df$Co[flip.df$Co=="ANGOLA" & df[p,"country_stand"]=="DEMOCRATIC REPUBLIC OF THE CONGO"]<-"DEMOCRATIC REPUBLIC OF THE CONGO"
				flip.df$Co[flip.df$Co=="OMAN" & df[p,"country_stand"]=="UNITED ARAB EMIRATES"]<-"UNITED ARAB EMIRATES"
				flip.df$Co[flip.df$Co=="PANAMA" & df[p,"country_stand"]=="COSTA RICA"]<-"COSTA RICA"
				flip.df$Co[flip.df$Co=="ARGENTINA" & df[p,"country_stand"]=="BRAZIL"]<-"BRAZIL"
				flip.df$Co[flip.df$Co=="PARAGUAY" & df[p,"country_stand"]=="ARGENTINA"]<-"ARGENTINA"
				flip.df$Co[flip.df$Co=="UNITED STATES" & df[p,"country_stand"]=="MEXICO"]<-"MEXICO"
				flip.df$Co[flip.df$Co=="SURINAME" & df[p,"country_stand"]=="GUYANA"]<-"GUYANA"
				flip.df$Co[flip.df$Co=="VENEZUELA" & df[p,"country_stand"]=="COLOMBIA"]<-"COLOMBIA"
				flip.df$Co[flip.df$Co=="CAYMAN ISLANDS" & df[p,"country_stand"]=="UNITED KINGDOM"]<-"UNITED KINGDOM"
				flip.df$Co[flip.df$Co=="PALAU" & df[p,"country_stand"]=="MICRONESIA"]<-"MICRONESIA"
				flip.df$Co[flip.df$Co=="COLOMBIA" & df[p,"country_stand"]=="PERU"]<-"PERU"
				flip.df$Co[flip.df$Co=="BRAZIL" & df[p,"country_stand"]=="BOLIVIA"]<-"BOLIVIA"
				flip.df$Co[flip.df$Co=="HAITI" & df[p,"country_stand"]=="UNITED STATES"]<-"UNITED STATES"
				flip.df$Co[flip.df$Co=="UNITED KINGDOM" & df[p,"country_stand"]=="SCOTLAND"]<-"SCOTLAND"
				flip.df$Co[flip.df$Co=="NETHERLANDS" & df[p,"country_stand"]=="BELGIUM"]<-"BELGIUM"
				flip.df$Co[flip.df$Co=="UNITED KINGDOM" & df[p,"country_stand"]=="ENGLAND"]<-"ENGLAND"
				flip.df$Co[flip.df$Co=="GERMANY" & df[p,"country_stand"]=="CZECH REPUBLIC"]<-"CZECH REPUBLIC"
				flip.df$Co[flip.df$Co=="TURKMENISTAN" & df[p,"country_stand"]=="IRAN"]<-"IRAN"
				flip.df$Co[flip.df$Co=="SOUTH AFRICA" & df[p,"country_stand"]=="BOTSWANA"]<-"BOTSWANA"
				flip.df$Co[flip.df$Co=="ETHIOPIA" & df[p,"country_stand"]=="SOMALIA"]<-"SOMALIA"
				flip.df$Co[flip.df$Co=="MADEIRA ISLANDS" & df[p,"country_stand"]=="PORTUGAL"]<-"PORTUGAL"
				flip.df$Co[flip.df$Co=="URUGUAY" & df[p,"country_stand"]=="ARGENTINA"]<-"ARGENTINA"
				flip.df$Co[flip.df$Co=="ARGENTINA" & df[p,"country_stand"]=="PARAGUAY"]<-"PARAGUAY"
				flip.df$Co[flip.df$Co=="AZORES" & df[p,"country_stand"]=="PORTUGAL"]<-"PORTUGAL"
				flip.df$Co[flip.df$Co=="VIRGIN ISLANDS, US" & df[p,"country_stand"]=="UNITED STATES VIRGIN ISLANDS"]<-"UNITED STATES VIRGIN ISLANDS"
				flip.df$Co[flip.df$Co=="GUAM" & df[p,"country_stand"]=="MICRONESIA"]<-"MICRONESIA"
				flip.df$Co[flip.df$Co=="CAMBODIA" & df[p,"country_stand"]=="VIETNAM"]<-"VIETNAM"
				flip.df$Co[flip.df$Co=="NORFOLK ISLAND" & df[p,"country_stand"]=="AUSTRALIA"]<-"AUSTRALIA"
				flip.df$Co[flip.df$Co=="SOUTH AFRICA" & df[p,"country_stand"]=="ZIMBABWE"]<-"ZIMBABWE"
				flip.df$Co[flip.df$Co=="SOUTH SUDAN" & df[p,"country_stand"]=="SUDAN"]<-"SUDAN"
				flip.df$Co[flip.df$Co=="BOLIVIA" & df[p,"country_stand"]=="PERU"]<-"PERU"
				flip.df$Co[flip.df$Co=="MALAYSIA" & df[p,"country_stand"]=="INDONESIA"]<-"INDONESIA"
				flip.df$Co[flip.df$Co=="INDONESIA" & df[p,"country_stand"]=="PAPUA NEW GUINEA"]<-"PAPUA NEW GUINEA"
				flip.df$Co[flip.df$Co=="EL SALVADOR" & df[p,"country_stand"]=="GUATEMALA"]<-"GUATEMALA"
				flip.df$Co[flip.df$Co=="BRUNEI" & df[p,"country_stand"]=="MALAYSIA"]<-"MALAYSIA"
				flip.df$Co[flip.df$Co=="MONTENEGRO" & df[p,"country_stand"]=="CROATIA"]<-"CROATIA"
				flip.df$Co[flip.df$Co=="BOLIVIA" & df[p,"country_stand"]=="ARGENTINA"]<-"ARGENTINA"
				flip.df$Co[flip.df$Co=="MYANMAR" & df[p,"country_stand"]=="THAILAND"]<-"THAILAND"
				flip.df$Co[flip.df$Co=="MEXICO" & df[p,"country_stand"]=="GUATEMALA"]<-"GUATEMALA"
				flip.df$Co[flip.df$Co=="PERU" & df[p,"country_stand"]=="ECUADOR"]<-"ECUADOR"
				flip.df$Co[flip.df$Co=="COLOMBIA" & df[p,"country_stand"]=="BRAZIL"]<-"BRAZIL"
				flip.df$Co[flip.df$Co=="BRAZIL" & df[p,"country_stand"]=="COLOMBIA"]<-"COLOMBIA"
				flip.df$Co[flip.df$Co=="LESOTHO" & df[p,"country_stand"]=="SOUTH AFRICA"]<-"SOUTH AFRICA"
				flip.df$Co[flip.df$Co=="ZAMBIA" & df[p,"country_stand"]=="ZIMBABWE"]<-"ZIMBABWE"
				flip.df$Co[flip.df$Co=="VIRGIN ISLANDS, BRITISH" & df[p,"country_stand"]=="BRITISH VIRGIN ISLANDS"]<-"BRITISH VIRGIN ISLANDS"
				flip.df$Co[flip.df$Co=="SAMOA" & df[p,"country_stand"]=="AMERICAN SAMOA"]<-"AMERICAN SAMOA"
				flip.df$Co[flip.df$Co=="REPUBLIC OF CONGO" & df[p,"country_stand"]=="CONGO"]<-"CONGO"
				flip.df$Co[flip.df$Co=="DEMOCRATIC REPUBLIC OF THE CONGO" & df[p,"country_stand"]=="CONGO"]<-"CONGO"
				#IF ANY OF THESE MATCH THE COUNTRY OF INTEREST...IF NOT, CHECK ORIGINAL WAS OVER LAND
				if(all(flip.df$Co!=df[p,"country_stand"],na.rm=TRUE)){
					df[p,c("decimal_longitude_corr","decimal_latitude_corr")]<-df[p,c("decimal_longitude","decimal_latitude")]
					df[p,"closest_country"]<-clos
				}	
				if(isTRUE(any(flip.df$Co==df[p,"country_stand"]))){
					df[p,c("decimal_longitude_corr","decimal_latitude_corr")]<-flip.df[!is.na(flip.df$Co) & flip.df$Co==df[p,"country_stand"],c("decimal_longitude","decimal_latitude")][1,]
					df[p,"closest_country"]<-df[p,"country_stand"]
					df[p,"coords_flipped"]<-"Yes"
				}
				df$Land[p]<-filterByLand(df[p,c("decimal_longitude_corr","decimal_latitude_corr")])
				if(isTRUE(df$Land[p]==TRUE & df$country_stand[p]==df$closest_country[p])){
					df$Keep[p]<-"Yes"
				}
				if(isTRUE(df$Land[p]==TRUE & is.na(df$closest_country[p]))){
					df$Keep[p]<-"Yes"
				}
				if(df$Keep[p]=="No"){
					if(isTRUE(df$Land[p]==TRUE)){
						print(paste(df$Taxon[p],df$decimal_longitude[p],df$decimal_latitude[p],"listed in", df$country_stand[p],"but matches",df$closest_country[p],sep=" "))
						zz<-nrow(problems)+1
						problems[zz,]<-c(df$Taxon[p],df$decimal_longitude[p],df$decimal_latitude[p],df$country_stand[p],df$closest_country[p],NA)
						problems$comb<-paste(problems$country_stand,problems$closest_country,sep="_")
						problems<-problems[!duplicated(problems$comb),]
					}
					if(isTRUE(df$Land[p]==FALSE)){
						#print(paste(df$Taxon[p],df$decimal_longitude[p],df$decimal_latitude[p],"listed in", df$country_stand[p],"but matches",df$closest_country[p],"and is not over land",sep=" "))
					}
				}
				rm(flip.df)
			}
		}
		write.csv(df,paste(file.path.mod,fls[x],"_bioclim_vars_TEST.csv",sep=""),row.names=FALSE)
		rm(clos)
		rm(df)
		gc()
	}
}	

write.csv(problems,"~/Desktop/ants_ecogeography/PROBLEMS_bioclim_vars_TEST.csv",row.names=FALSE)



bc<-paste("bio",1:19,sep="_")
bc<-c(bc,"NPP","PET","ELEV","BIOME","REALM")
my.path<-"~/Desktop/ants_ecogeography/new_enm/ascii_layers_for_maxent/"
my.filenames<-paste(my.path,bc,".asc",sep="")

file.path.mod<-"~/Desktop/ants_ecogeography/antweb_bioclim_vars_06sep2017/individual_files_w_addl_data_countrymods/"
fls<-list.files(path=file.path.mod)
fls<-fls[grep("_bioclim_vars_TEST.csv",fls)]
fls<-gsub("_bioclim_vars_TEST.csv","",fls)

uniques<-data.frame(matrix(nrow=0,ncol=51))
colnames(uniques)<-c("url","catalogNumber","family","subfamily","genus","specificEpithet","scientific_name","typeStatus","stateProvince","country","dateIdentified","habitat","minimumElevationInMeters","geojson.type","decimal_latitude","decimal_longitude","bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8","bio9","bio10","bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18","bio19","Taxon","Biome","AnnTotPrecip","AvgAnnRH","Elev","Evapotrans","NPP","PotEvapotrans","country_mod","country_stand","closest_country","decimal_longitude_corr","decimal_latitude_corr","coords_flipped","Land","Keep")
#zz<-0

drops<-data.frame(matrix(nrow=0,ncol=51))
colnames(drops)<-c("url","catalogNumber","family","subfamily","genus","specificEpithet","scientific_name","typeStatus","stateProvince","country","dateIdentified","habitat","minimumElevationInMeters","geojson.type","decimal_latitude","decimal_longitude","bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8","bio9","bio10","bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18","bio19","Taxon","Biome","AnnTotPrecip","AvgAnnRH","Elev","Evapotrans","NPP","PotEvapotrans","country_mod","country_stand","closest_country","decimal_longitude_corr","decimal_latitude_corr","coords_flipped","Land","Keep")

#get drops
for(x in 1:length(fls)){
	df<-read.csv(file=paste(file.path.mod,fls[x],"_bioclim_vars_TEST.csv",sep=""),stringsAsFactors=FALSE)
	df<-df[df$Keep=="No",]
	if(nrow(df)>0){
		#remove duplicated lat/longs and add to uniques
		drops<-rbind(drops,df[!duplicated(df[,c("decimal_longitude_corr","decimal_latitude_corr")]),])
		#remove duplicated lat/longs from uniques
		#uniques<-uniques[!duplicated(uniques[,c("decimal_longitude_corr","decimal_latitude_corr")]),]
	}
}

write.csv(drops,file="~/Desktop/ants_ecogeography/drops_coords_per_species.csv",row.names=FALSE)


#get unique lat/longs from all all retained
for(x in 1:length(fls)){
	df<-read.csv(file=paste(file.path.mod,fls[x],"_bioclim_vars_TEST.csv",sep=""),stringsAsFactors=FALSE)
	df<-df[df$Keep=="Yes",]
	if(nrow(df)>0){
		#remove duplicated lat/longs and add to uniques
		uniques<-rbind(uniques,df[!duplicated(df[,c("decimal_longitude_corr","decimal_latitude_corr")]),])
		#remove duplicated lat/longs from uniques
		#uniques<-uniques[!duplicated(uniques[,c("decimal_longitude_corr","decimal_latitude_corr")]),]
	}
}

write.csv(uniques,file="~/Desktop/ants_ecogeography/unique_coords_per_species.csv",row.names=FALSE)


#essentially taken from here: http://stackoverflow.com/questions/27338512/color-countries-on-world-map-based-on-iso3-codes-in-r-using-ggplot
library(maptools)
library(mapproj)
library(rgeos)
library(rgdal)
library(ggplot2)
library(jsonlite)
library(RCurl)
library(raster)
library(RColorBrewer)

uniques<-read.csv(file="~/Desktop/ants_ecogeography/unique_coords_per_species.csv",stringsAsFactors=FALSE)
cols<-colnames(uniques)
keepers<-cols[c(1:16,36,44:51)]
uniques.trim<-uniques[,keepers]

#from http://stackoverflow.com/questions/14334970/convert-latitude-and-longitude-coordinates-to-country-name-in-r
biomes<-readOGR(dsn="/Users/matthewnelsen/Desktop/ants_ecogeography/official",layer="wwf_terr_ecos")
#98 and 99 represent water or other things so replace with NA
biomes$BIOME[biomes$BIOME==98]<-NA
biomes$BIOME[biomes$BIOME==99]<-NA
unique(biomes$BIOME)

for(x in 0:3){
	if(x<3){
		temp<-extract(x=biomes,y=cbind(uniques.trim$decimal_longitude_corr[((50000*x)+1):(50000*(x+1))],uniques.trim$decimal_latitude_corr[((50000*x)+1):(50000*(x+1))]),method="simple")
		uniques.trim[((50000*x)+1):(50000*(x+1)),"BIOME"]<-temp[,"BIOME"]
		uniques.trim[((50000*x)+1):(50000*(x+1)),"REALM"]<-temp[,"REALM"]		
	}
	if(x==3){
		temp<-extract(x=biomes,y=cbind(uniques.trim$decimal_longitude_corr[((50000*x)+1):nrow(uniques.trim)],uniques.trim$decimal_latitude_corr[((50000*x)+1):nrow(uniques.trim)]),method="simple")
		uniques.trim[((50000*x)+1):nrow(uniques.trim),"BIOME"]<-temp[,"BIOME"]
		uniques.trim[((50000*x)+1):nrow(uniques.trim),"REALM"]<-temp[,"REALM"]		
	}
}

npp.file<-"~/Desktop/ants_ecogeography/new_enm/MOD17A3_Science_NPP_mean_00_15.tif"
npp<-raster(npp.file,values=TRUE)
npp[npp==65535]<-NA
npp

for(x in 0:3){
	if(x<3){
		uniques.trim[((50000*x)+1):(50000*(x+1)),"NPP"]<-extract(x=npp,y=cbind(uniques.trim$decimal_longitude_corr[((50000*x)+1):(50000*(x+1))],uniques.trim$decimal_latitude_corr[((50000*x)+1):(50000*(x+1))]),method="simple")
	}
	if(x==3){
		uniques.trim[((50000*x)+1):nrow(uniques.trim),"NPP"]<-extract(x=npp,y=cbind(uniques.trim$decimal_longitude_corr[((50000*x)+1):nrow(uniques.trim)],uniques.trim$decimal_latitude_corr[((50000*x)+1):nrow(uniques.trim)]),method="simple")
	}
}

#from https://gis.stackexchange.com/questions/132403/how-to-read-adf-files-into-r
dpath<-"~/Desktop/ants_ecogeography/new_enm/Global PET and Aridity Index/PET_he_annual/pet_he_yr"
x<-new("GDALReadOnlyDataset",dpath)
getDriver(x)
getDriverLongName(getDriver(x))
xx<-asSGDF_GROD(x)
pet <- raster(xx)
pet

for(x in 0:3){
	if(x<3){
		uniques.trim[((50000*x)+1):(50000*(x+1)),"PET"]<-extract(x=pet,y=cbind(uniques.trim$decimal_longitude_corr[((50000*x)+1):(50000*(x+1))],uniques.trim$decimal_latitude_corr[((50000*x)+1):(50000*(x+1))]),method="simple")
	}
	if(x==3){
		uniques.trim[((50000*x)+1):nrow(uniques.trim),"PET"]<-extract(x=pet,y=cbind(uniques.trim$decimal_longitude_corr[((50000*x)+1):nrow(uniques.trim)],uniques.trim$decimal_latitude_corr[((50000*x)+1):nrow(uniques.trim)]),method="simple")
	}
}

#GET ELEVATION AND ADD
#from https://gis.stackexchange.com/questions/132403/how-to-read-adf-files-into-r
dpath2<-"~/Desktop/ants_ecogeography/new_enm/gis_layers/alt_30s_esri/alt/alt"
xyz<-new("GDALReadOnlyDataset",dpath2)
getDriver(xyz)
getDriverLongName(getDriver(xyz))
xxyz<-asSGDF_GROD(xyz)
elev <- raster(xxyz)
for(x in 0:3){
	if(x<3){
		uniques.trim[((50000*x)+1):(50000*(x+1)),"ELEV"]<-extract(x=elev,y=cbind(uniques.trim$decimal_longitude_corr[((50000*x)+1):(50000*(x+1))],uniques.trim$decimal_latitude_corr[((50000*x)+1):(50000*(x+1))]),method="simple")
	}
	if(x==3){
		uniques.trim[((50000*x)+1):nrow(uniques.trim),"ELEV"]<-extract(x=elev,y=cbind(uniques.trim$decimal_longitude_corr[((50000*x)+1):nrow(uniques.trim)],uniques.trim$decimal_latitude_corr[((50000*x)+1):nrow(uniques.trim)]),method="simple")
	}
}

ls<-getData(name="worldclim",var="bio", res=2.5, download=TRUE)
uniques.trim[,paste("bio",1:19,"orig",sep="_")]<-NA
for(x in 0:3){
	if(x<3){
		uniques.trim[((50000*x)+1):(50000*(x+1)),paste("bio",1:19,"orig",sep="_")]<-extract(x=ls,y=cbind(uniques.trim$decimal_longitude_corr[((50000*x)+1):(50000*(x+1))],uniques.trim$decimal_latitude_corr[((50000*x)+1):(50000*(x+1))]),method="simple")
	}
	if(x==3){
		uniques.trim[((50000*x)+1):nrow(uniques.trim),paste("bio",1:19,"orig",sep="_")]<-extract(x=ls,y=cbind(uniques.trim$decimal_longitude_corr[((50000*x)+1):nrow(uniques.trim)],uniques.trim$decimal_latitude_corr[((50000*x)+1):nrow(uniques.trim)]),method="simple")
	}
}

sum(is.na(uniques.trim$PET))
sum(is.na(uniques.trim$ELEV))
sum(is.na(uniques.trim$BIOME))
sum(is.na(uniques.trim$bio_1_orig))
sum(rowSums(is.na(uniques.trim[,c("BIOME","REALM","NPP","PET","ELEV",paste("bio",1:19,"orig",sep="_"))]))>0)

length(unique(uniques.trim$Taxon))
length(unique(uniques.trim$Taxon[rowSums(is.na(uniques.trim[,c("BIOME","REALM","NPP","PET","ELEV",paste("bio",1:19,"orig",sep="_"))]))>0]))

write.csv(uniques.trim,file="~/Desktop/ants_ecogeography/unique_coords_per_species_w_dat_24May2018.csv",row.names=FALSE)

tectonics<-readOGR(dsn="/Users/matthewnelsen/Desktop/tectonicplates-master",layer="PB2002_plates_MODIFIED_Third")
unique(tectonics$PlateName)

uniques.trim<-read.csv(file="~/Desktop/ants_ecogeography/unique_coords_per_species_w_dat_24May2018.csv",stringsAsFactors=FALSE)

for(x in 0:3){
	if(x<3){
		temp<-extract(x=tectonics,y=cbind(uniques.trim$decimal_longitude_corr[((50000*x)+1):(50000*(x+1))],uniques.trim$decimal_latitude_corr[((50000*x)+1):(50000*(x+1))]),method="simple")
		uniques.trim[((50000*x)+1):(50000*(x+1)),"Code"]<-temp[,"Code"]
		uniques.trim[((50000*x)+1):(50000*(x+1)),"PlateName"]<-temp[,"PlateName"]		
	}
	if(x==3){
		temp<-extract(x=tectonics,y=cbind(uniques.trim$decimal_longitude_corr[((50000*x)+1):nrow(uniques.trim)],uniques.trim$decimal_latitude_corr[((50000*x)+1):nrow(uniques.trim)]),method="simple")
		uniques.trim[((50000*x)+1):nrow(uniques.trim),"Code"]<-temp[,"Code"]
		uniques.trim[((50000*x)+1):nrow(uniques.trim),"PlateName"]<-temp[,"PlateName"]		
	}
}
uniques.trim$Code<-as.character(uniques.trim$Code)
uniques.trim$PlateName<-as.character(uniques.trim$PlateName)
#uniques.trim$Code[uniques.trim$Code=="NA"]<-"NO"

write.csv(uniques.trim,file="~/Desktop/ants_ecogeography/unique_coords_per_species_w_dat_plates_24May2018.csv",row.names=FALSE)

#above was conducted as used_this8.  Now am adding in soil variables
require(raster)
t_gravel_file <- '/Volumes/welwitschia/hwsd_grds/t_gravel.grd'
t_sand_file <- '/Volumes/welwitschia/hwsd_grds/t_sand.grd'
t_silt_file <- '/Volumes/welwitschia/hwsd_grds/t_silt.grd'
t_clay_file <- '/Volumes/welwitschia/hwsd_grds/t_clay.grd'
t_bulk_density_file <- '/Volumes/welwitschia/hwsd_grds/t_bulk_density.grd'
t_oc_file <- '/Volumes/welwitschia/hwsd_grds/t_oc.grd'
t_ph_h20_file <- '/Volumes/welwitschia/hwsd_grds/t_ph_h20.grd'
t_cec_clay_file <- '/Volumes/welwitschia/hwsd_grds/t_cec_clay.grd'
t_cec_soil_file <- '/Volumes/welwitschia/hwsd_grds/t_cec_soil.grd'
t_bs_file <- '/Volumes/welwitschia/hwsd_grds/t_bs.grd'
t_teb_file <- '/Volumes/welwitschia/hwsd_grds/t_teb.grd'
t_caco3_file <- '/Volumes/welwitschia/hwsd_grds/t_caco3.grd'
t_caso4_file <- '/Volumes/welwitschia/hwsd_grds/t_caso4.grd'
t_esp_file <- '/Volumes/welwitschia/hwsd_grds/t_esp.grd'
t_ece_file <- '/Volumes/welwitschia/hwsd_grds/t_ece.grd'

T_GRAVEL <- raster(t_gravel_file)
T_SAND <- raster(t_sand_file)
T_SILT <- raster(t_silt_file)
T_CLAY <- raster(t_clay_file)
T_BULK_DENSITY <- raster(t_bulk_density_file)
T_OC <- raster(t_oc_file)
T_PH_H2O <- raster(t_ph_h20_file)
T_CEC_CLAY <- raster(t_cec_clay_file)
T_CEC_SOIL <- raster(t_cec_soil_file)
T_BS <- raster(t_bs_file)
T_TEB <- raster(t_teb_file)
T_CACO3 <- raster(t_caco3_file)
T_CASO4 <- raster(t_caso4_file)
T_ESP <- raster(t_esp_file)
T_ECE <- raster(t_ece_file)

new.output<-read.csv(file="~/Desktop/ants_ecogeography/unique_coords_per_species_w_dat_plates_24May2018.csv",stringsAsFactors=FALSE)
new.output$Taxon[new.output$Taxon=="Gesomyrmex_kh01"]<-"Gesomyrmex_CASENT0106178"
new.output$Taxon[new.output$Taxon=="Melophorus_au01"]<-"Melophorus_CASENT0106148"
new.output$Taxon[new.output$Taxon=="Protanilla_jp01"]<-"Protanilla_CASENT0007002"
new.output$Taxon[new.output$Taxon=="Leptanilla_gr01"]<-"Leptanilla_CASENT0106505"

cts<-c("T_GRAVEL","T_SAND","T_SILT","T_CLAY","T_BULK_DENSITY","T_OC","T_PH_H2O","T_CEC_CLAY","T_CEC_SOIL","T_BS","T_TEB","T_CACO3","T_CASO4","T_ESP","T_ECE")

new.output[,cts]<-NA

for(x in 0:3){
	if(x<3){
		new.output[((50000*x)+1):(50000*(x+1)),"T_GRAVEL"]<-extract(x=T_GRAVEL,y=cbind(new.output$decimal_longitude_corr[((50000*x)+1):(50000*(x+1))], new.output$decimal_latitude_corr[((50000*x)+1):(50000*(x+1))]),method="simple")
		new.output[((50000*x)+1):(50000*(x+1)),"T_SAND"]<-extract(x=T_SAND,y=cbind(new.output$decimal_longitude_corr[((50000*x)+1):(50000*(x+1))], new.output$decimal_latitude_corr[((50000*x)+1):(50000*(x+1))]),method="simple")
		new.output[((50000*x)+1):(50000*(x+1)),"T_SILT"]<-extract(x=T_SILT,y=cbind(new.output$decimal_longitude_corr[((50000*x)+1):(50000*(x+1))], new.output$decimal_latitude_corr[((50000*x)+1):(50000*(x+1))]),method="simple")
		new.output[((50000*x)+1):(50000*(x+1)),"T_CLAY"]<-extract(x=T_CLAY,y=cbind(new.output$decimal_longitude_corr[((50000*x)+1):(50000*(x+1))], new.output$decimal_latitude_corr[((50000*x)+1):(50000*(x+1))]),method="simple")
		new.output[((50000*x)+1):(50000*(x+1)),"T_BULK_DENSITY"]<-extract(x=T_BULK_DENSITY,y=cbind(new.output$decimal_longitude_corr[((50000*x)+1):(50000*(x+1))], new.output$decimal_latitude_corr[((50000*x)+1):(50000*(x+1))]),method="simple")
		new.output[((50000*x)+1):(50000*(x+1)),"T_OC"]<-extract(x=T_OC,y=cbind(new.output$decimal_longitude_corr[((50000*x)+1):(50000*(x+1))], new.output$decimal_latitude_corr[((50000*x)+1):(50000*(x+1))]),method="simple")
		new.output[((50000*x)+1):(50000*(x+1)),"T_PH_H2O"]<-extract(x=T_PH_H2O,y=cbind(new.output$decimal_longitude_corr[((50000*x)+1):(50000*(x+1))], new.output$decimal_latitude_corr[((50000*x)+1):(50000*(x+1))]),method="simple")
		new.output[((50000*x)+1):(50000*(x+1)),"T_CEC_CLAY"]<-extract(x=T_CEC_CLAY,y=cbind(new.output$decimal_longitude_corr[((50000*x)+1):(50000*(x+1))], new.output$decimal_latitude_corr[((50000*x)+1):(50000*(x+1))]),method="simple")
		new.output[((50000*x)+1):(50000*(x+1)),"T_CEC_SOIL"]<-extract(x=T_CEC_SOIL,y=cbind(new.output$decimal_longitude_corr[((50000*x)+1):(50000*(x+1))], new.output$decimal_latitude_corr[((50000*x)+1):(50000*(x+1))]),method="simple")
		new.output[((50000*x)+1):(50000*(x+1)),"T_BS"]<-extract(x=T_BS,y=cbind(new.output$decimal_longitude_corr[((50000*x)+1):(50000*(x+1))], new.output$decimal_latitude_corr[((50000*x)+1):(50000*(x+1))]),method="simple")
		new.output[((50000*x)+1):(50000*(x+1)),"T_TEB"]<-extract(x=T_TEB,y=cbind(new.output$decimal_longitude_corr[((50000*x)+1):(50000*(x+1))], new.output$decimal_latitude_corr[((50000*x)+1):(50000*(x+1))]),method="simple")
		new.output[((50000*x)+1):(50000*(x+1)),"T_CACO3"]<-extract(x=T_CACO3,y=cbind(new.output$decimal_longitude_corr[((50000*x)+1):(50000*(x+1))], new.output$decimal_latitude_corr[((50000*x)+1):(50000*(x+1))]),method="simple")
		new.output[((50000*x)+1):(50000*(x+1)),"T_CASO4"]<-extract(x=T_CASO4,y=cbind(new.output$decimal_longitude_corr[((50000*x)+1):(50000*(x+1))], new.output$decimal_latitude_corr[((50000*x)+1):(50000*(x+1))]),method="simple")
		new.output[((50000*x)+1):(50000*(x+1)),"T_ESP"]<-extract(x= T_ESP,y=cbind(new.output$decimal_longitude_corr[((50000*x)+1):(50000*(x+1))], new.output$decimal_latitude_corr[((50000*x)+1):(50000*(x+1))]),method="simple")
		new.output[((50000*x)+1):(50000*(x+1)),"T_ECE"]<-extract(x=T_ECE,y=cbind(new.output$decimal_longitude_corr[((50000*x)+1):(50000*(x+1))], new.output$decimal_latitude_corr[((50000*x)+1):(50000*(x+1))]),method="simple")
	}
	if(x==3){
		new.output[((50000*x)+1):nrow(new.output),"T_GRAVEL"]<-extract(x=T_GRAVEL,y=cbind(new.output$decimal_longitude_corr[((50000*x)+1):nrow(new.output)], new.output$decimal_latitude_corr[((50000*x)+1):nrow(new.output)]),method="simple")
		new.output[((50000*x)+1):nrow(new.output),"T_SAND"]<-extract(x=T_SAND,y=cbind(new.output$decimal_longitude_corr[((50000*x)+1):nrow(new.output)], new.output$decimal_latitude_corr[((50000*x)+1):nrow(new.output)]),method="simple")
		new.output[((50000*x)+1):nrow(new.output),"T_SILT"]<-extract(x=T_SILT,y=cbind(new.output$decimal_longitude_corr[((50000*x)+1):nrow(new.output)], new.output$decimal_latitude_corr[((50000*x)+1):nrow(new.output)]),method="simple")
		new.output[((50000*x)+1):nrow(new.output),"T_CLAY"]<-extract(x=T_CLAY,y=cbind(new.output$decimal_longitude_corr[((50000*x)+1):nrow(new.output)], new.output$decimal_latitude_corr[((50000*x)+1):nrow(new.output)]),method="simple")
		new.output[((50000*x)+1):nrow(new.output),"T_BULK_DENSITY"]<-extract(x=T_BULK_DENSITY,y=cbind(new.output$decimal_longitude_corr[((50000*x)+1):nrow(new.output)], new.output$decimal_latitude_corr[((50000*x)+1):nrow(new.output)]),method="simple")
		new.output[((50000*x)+1):nrow(new.output),"T_OC"]<-extract(x=T_OC,y=cbind(new.output$decimal_longitude_corr[((50000*x)+1):nrow(new.output)], new.output$decimal_latitude_corr[((50000*x)+1):nrow(new.output)]),method="simple")
		new.output[((50000*x)+1):nrow(new.output),"T_PH_H2O"]<-extract(x=T_PH_H2O,y=cbind(new.output$decimal_longitude_corr[((50000*x)+1):nrow(new.output)], new.output$decimal_latitude_corr[((50000*x)+1):nrow(new.output)]),method="simple")
		new.output[((50000*x)+1):nrow(new.output),"T_CEC_CLAY"]<-extract(x=T_CEC_CLAY,y=cbind(new.output$decimal_longitude_corr[((50000*x)+1):nrow(new.output)], new.output$decimal_latitude_corr[((50000*x)+1):nrow(new.output)]),method="simple")
		new.output[((50000*x)+1):nrow(new.output),"T_CEC_SOIL"]<-extract(x=T_CEC_SOIL,y=cbind(new.output$decimal_longitude_corr[((50000*x)+1):nrow(new.output)], new.output$decimal_latitude_corr[((50000*x)+1):nrow(new.output)]),method="simple")
		new.output[((50000*x)+1):nrow(new.output),"T_BS"]<-extract(x=T_BS,y=cbind(new.output$decimal_longitude_corr[((50000*x)+1):nrow(new.output)], new.output$decimal_latitude_corr[((50000*x)+1):nrow(new.output)]),method="simple")
		new.output[((50000*x)+1):nrow(new.output),"T_TEB"]<-extract(x=T_TEB,y=cbind(new.output$decimal_longitude_corr[((50000*x)+1):nrow(new.output)], new.output$decimal_latitude_corr[((50000*x)+1):nrow(new.output)]),method="simple")
		new.output[((50000*x)+1):nrow(new.output),"T_CACO3"]<-extract(x=T_CACO3,y=cbind(new.output$decimal_longitude_corr[((50000*x)+1):nrow(new.output)], new.output$decimal_latitude_corr[((50000*x)+1):nrow(new.output)]),method="simple")
		new.output[((50000*x)+1):nrow(new.output),"T_CASO4"]<-extract(x=T_CASO4,y=cbind(new.output$decimal_longitude_corr[((50000*x)+1):nrow(new.output)], new.output$decimal_latitude_corr[((50000*x)+1):nrow(new.output)]),method="simple")
		new.output[((50000*x)+1):nrow(new.output),"T_ESP"]<-extract(x=T_ESP,y=cbind(new.output$decimal_longitude_corr[((50000*x)+1):nrow(new.output)], new.output$decimal_latitude_corr[((50000*x)+1):nrow(new.output)]),method="simple")
		new.output[((50000*x)+1):nrow(new.output),"T_ECE"]<-extract(x=T_ECE,y=cbind(new.output$decimal_longitude_corr[((50000*x)+1):nrow(new.output)], new.output$decimal_latitude_corr[((50000*x)+1):nrow(new.output)]),method="simple")	}
}
write.csv(new.output,file="~/Desktop/ants_ecogeography/unique_coords_per_species_w_dat_plates_24May2018_soils_added.csv",row.names=FALSE)



new.output<-read.csv(file="~/Desktop/ants_ecogeography/unique_coords_per_species_w_dat_plates_24May2018_soils_added.csv",stringsAsFactors=FALSE)

#drop T_CEC_CLAY as it refers to clay portion and many lack clay-rich soil and are NA
soil.cts<-c("T_GRAVEL","T_SAND","T_SILT","T_CLAY","T_BULK_DENSITY","T_OC","T_PH_H2O","T_CEC_SOIL","T_BS","T_TEB","T_CACO3","T_CASO4","T_ESP","T_ECE")
env.cts<-c("BIOME","NPP","PET","ELEV",paste("bio",1:19,"orig",sep="_"))
cts<-c(env.cts,soil.cts)

#REMOVE ANY W NA's
new.output.nona<-new.output[rowSums(is.na(new.output[,cts]))==0,]
nrow(new.output.nona)
new.output.nona.env.only<-new.output[rowSums(is.na(new.output[,env.cts]))==0,]
nrow(new.output.nona.env.only)
new.output.nona.soil.only<-new.output[rowSums(is.na(new.output[,soil.cts]))==0,]
nrow(new.output.nona.soil.only)


intro<-read.csv(file="~/Desktop/ants_ecogeography/antweb.introduced.28dec2016.clean.csv",stringsAsFactors=FALSE)
new.output.no.intr<-new.output.nona[!new.output.nona$Taxon %in% intro$Taxon,]

new.tab.cols<-sort(c("AN",unique(new.output.no.intr$Code)))
new.tab<-as.data.frame(matrix(nrow=length(unique(new.output.no.intr$Taxon)),ncol=length(new.tab.cols)+1))
rownames(new.tab)<-unique(new.output.no.intr$Taxon)
colnames(new.tab)<-c(new.tab.cols,"TOTAL")
new.tab[,c(1:ncol(new.tab))]<-0
for(x in 1:nrow(new.tab)){
	plate.tab<-NA
	plate.tab<-table(new.output.no.intr$Code[new.output.no.intr$Taxon %in% rownames(new.tab)[x]])
	for(z in 1:length(plate.tab)){
		new.tab[x,colnames(new.tab) %in% names(plate.tab)[z]]<-plate.tab[z][[1]]
	}
	new.tab$TOTAL[x]<-rowSums(new.tab[x,])	
}

new.tab.prop<-new.tab[,c(1:(ncol(new.tab)-1))]
new.tab.prop[,1:ncol(new.tab.prop)]<-0
for(x in 1:nrow(new.tab.prop)){
	for(z in 1:ncol(new.tab.prop)){
		new.tab.prop[x,z]<-100*(new.tab[x,z]/new.tab[x,"TOTAL"])
	}
}
new.tab.prop.min<-new.tab.prop
new.tab.prop.min[new.tab.prop.min<15]<-0
new.tab.prop.min[new.tab.prop.min>=15]<-1

max(rowSums(new.tab.prop.min))
min(rowSums(new.tab.prop.min))
#some have more than three regions

#just keep taxa in tree
require(phytools)
ants.tr<-read.tree(file="~/Desktop/14july2016_antcat_aligns_fullalign_einsi_nogappy_wadds_heavilymodified_correctedmacse/Hymenoptera/drop_bads_again/drop_bads_again_partfinder/proper_raxml_ml500searches_Pristo_Out/formicidae_ml_treepl_147/bamm_analyses_185_10nov2016/ML_TREE_treepl_185_outsdroppeds.ladderized.increasing.tre")

keepers<-new.tab.prop.min[rownames(new.tab.prop.min) %in% ants.tr$tip.label,]
nrow(keepers)
max(rowSums(keepers))
min(rowSums(keepers))

new.tab.prop[rownames(new.tab.prop) %in% rownames(keepers[rowSums(keepers)>3,]),]
#save this presence/absence for biogeography analysis


#now prune new.output.no.intr to only include taxa in tree and w greater than or equal to 15% in region
new.output.no.intr.prekeeps<-new.output.no.intr[new.output.no.intr$Taxon %in% rownames(keepers),]
new.output.no.intr.prekeeps$KeepAgain<-"No"

for(x in 1:nrow(keepers)){
	new.output.no.intr.prekeeps$KeepAgain[new.output.no.intr.prekeeps$Taxon %in% rownames(keepers)[x] & new.output.no.intr.prekeeps$Code %in% colnames(keepers)[keepers[x,]==1]]<-"Yes"
}

new.output.no.intr.keeps<-new.output.no.intr.prekeeps[new.output.no.intr.prekeeps$KeepAgain=="Yes",]


new.tabs<-as.data.frame.matrix(table(new.output.no.intr.keeps$Taxon,new.output.no.intr.keeps$Code))
new.tabs$AN<-0
new.tabs<-new.tabs[,order(colnames(new.tabs))]

min(rowSums(new.tabs))


new.tabs.bin<-new.tabs
new.tabs.bin[new.tabs.bin<1]<-0
new.tabs.bin[new.tabs.bin>=1]<-1
max(rowSums(new.tabs.bin))
min(rowSums(new.tabs.bin))


head(keepers)
head(new.tabs.bin)
identical(keepers,new.tabs.bin)

drops<-ants.tr$tip.label[!ants.tr$tip.label %in% rownames(new.tabs.bin)]
ants.tr.clean<-drop.tip(ants.tr,drops)


write.table(new.tabs.bin,file="~/Desktop/hab.data.plates.min.thresh29may2018.wsoils.nona.txt",quote=FALSE,sep="\t",col.names=NA)
write.tree(ants.tr.clean,file="~/Desktop/hab.data.plates.min.thresh.29may2018.wsoils.nona.tre")


write.csv(new.output.no.intr.keeps,file="~/Desktop/ants_ecogeography/unique_coords_per_species_cleaned_29May2018.wsoils.nona.csv")


sum.df<-read.csv(file="~/Desktop/ants_ecogeography/unique_coords_per_species_cleaned_29May2018.wsoils.nona.csv",stringsAsFactors=FALSE)
meds<-aggregate(sum.df[c("NPP","PET","ELEV",paste("bio",1:19,"orig",sep="_"),"T_GRAVEL","T_SAND","T_SILT","T_CLAY","T_BULK_DENSITY","T_OC","T_PH_H2O","T_CEC_SOIL","T_BS","T_TEB","T_CACO3","T_CASO4","T_ESP","T_ECE")], list(sum.df$Taxon),median,na.rm=TRUE)
colnames(meds)<-c("Taxon","NPP_Median","PET_Median","ELEV_Median",paste("bio",1:19,"orig_Median",sep="_"),paste(c("T_GRAVEL","T_SAND","T_SILT","T_CLAY","T_BULK_DENSITY","T_OC","T_PH_H2O","T_CEC_SOIL","T_BS","T_TEB","T_CACO3","T_CASO4","T_ESP","T_ECE"),"_Median",sep=""))
#mean
means<-aggregate(sum.df[c("NPP","PET","ELEV",paste("bio",1:19,"orig",sep="_"),"T_GRAVEL","T_SAND","T_SILT","T_CLAY","T_BULK_DENSITY","T_OC","T_PH_H2O","T_CEC_SOIL","T_BS","T_TEB","T_CACO3","T_CASO4","T_ESP","T_ECE")], list(sum.df$Taxon),mean,na.rm=TRUE)
colnames(means)<-c("Taxon","NPP_Mean","PET_Mean","ELEV_Mean",paste("bio",1:19,"orig_Mean",sep="_"),paste(c("T_GRAVEL","T_SAND","T_SILT","T_CLAY","T_BULK_DENSITY","T_OC","T_PH_H2O","T_CEC_SOIL","T_BS","T_TEB","T_CACO3","T_CASO4","T_ESP","T_ECE"),"_Mean",sep=""))
#sd
sd.samps<-aggregate(sum.df[c("NPP","PET","ELEV",paste("bio",1:19,"orig",sep="_"),"T_GRAVEL","T_SAND","T_SILT","T_CLAY","T_BULK_DENSITY","T_OC","T_PH_H2O","T_CEC_SOIL","T_BS","T_TEB","T_CACO3","T_CASO4","T_ESP","T_ECE")], list(sum.df$Taxon),sd,na.rm=TRUE)
colnames(sd.samps)<-c("Taxon","NPP_SD","PET_SD","ELEV_SD",paste("bio",1:19,"orig_SD",sep="_"),paste(c("T_GRAVEL","T_SAND","T_SILT","T_CLAY","T_BULK_DENSITY","T_OC","T_PH_H2O","T_CEC_SOIL","T_BS","T_TEB","T_CACO3","T_CASO4","T_ESP","T_ECE"),"_SD",sep=""))

#elevation range
max.meas<-aggregate(sum.df[c("NPP","PET","ELEV",paste("bio",1:19,"orig",sep="_"),"T_GRAVEL","T_SAND","T_SILT","T_CLAY","T_BULK_DENSITY","T_OC","T_PH_H2O","T_CEC_SOIL","T_BS","T_TEB","T_CACO3","T_CASO4","T_ESP","T_ECE")], list(sum.df$Taxon),max,na.rm=TRUE)
colnames(max.meas)<-c("Taxon","NPP_Max","PET_Max","ELEV_Max",paste("bio",1:19,"orig_Max",sep="_"),paste(c("T_GRAVEL","T_SAND","T_SILT","T_CLAY","T_BULK_DENSITY","T_OC","T_PH_H2O","T_CEC_SOIL","T_BS","T_TEB","T_CACO3","T_CASO4","T_ESP","T_ECE"),"_Max",sep=""))
max.meas[max.meas[,"NPP_Max"]==-Inf,"NPP_Max"]<-NA
max.meas[max.meas[,"PET_Max"]==-Inf,"PET_Max"]<-NA
max.meas[max.meas[,"ELEV_Max"]==-Inf,"ELEV_Max"]<-NA
min.meas<-aggregate(sum.df[c("NPP","PET","ELEV",paste("bio",1:19,"orig",sep="_"),"T_GRAVEL","T_SAND","T_SILT","T_CLAY","T_BULK_DENSITY","T_OC","T_PH_H2O","T_CEC_SOIL","T_BS","T_TEB","T_CACO3","T_CASO4","T_ESP","T_ECE")], list(sum.df$Taxon),min,na.rm=TRUE)
colnames(min.meas)<-c("Taxon","NPP_Min","PET_Min","ELEV_Min",paste("bio",1:19,"orig_Min",sep="_"),paste(c("T_GRAVEL","T_SAND","T_SILT","T_CLAY","T_BULK_DENSITY","T_OC","T_PH_H2O","T_CEC_SOIL","T_BS","T_TEB","T_CACO3","T_CASO4","T_ESP","T_ECE"),"_Min",sep=""))
min.meas[min.meas[,"NPP_Min"]==-Inf,"NPP_Min"]<-NA
min.meas[min.meas[,"PET_Min"]==-Inf,"PET_Min"]<-NA
min.meas[min.meas[,"ELEV_Min"]==-Inf,"ELEV_Min"]<-NA
ranges<-max.meas[,c(2:ncol(max.meas))]-min.meas[,c(2:ncol(min.meas))]
ranges<-cbind(max.meas$Taxon,ranges)
colnames(ranges)<-c("Taxon","NPP_Range","PET_Range","ELEV_Range",paste("bio",1:19,"orig_Range",sep="_"),paste(c("T_GRAVEL","T_SAND","T_SILT","T_CLAY","T_BULK_DENSITY","T_OC","T_PH_H2O","T_CEC_SOIL","T_BS","T_TEB","T_CACO3","T_CASO4","T_ESP","T_ECE"),"_Range",sep=""))

#counts those that are not NA
counts<-aggregate(sum.df[c("NPP","PET","ELEV",paste("bio",1,"orig",sep="_"),"BIOME","Code","T_GRAVEL","T_SAND","T_SILT","T_CLAY","T_BULK_DENSITY","T_OC","T_PH_H2O","T_CEC_SOIL","T_BS","T_TEB","T_CACO3","T_CASO4","T_ESP","T_ECE")], list(sum.df$Taxon),function(x) sum(!is.na(x)))
colnames(counts)<-c("Taxon","NPP_Count","PET_Count","ELEV_Count","bio_orig_Count","BIOME_Count","Code_Count",paste(c("T_GRAVEL","T_SAND","T_SILT","T_CLAY","T_BULK_DENSITY","T_OC","T_PH_H2O","T_CEC_SOIL","T_BS","T_TEB","T_CACO3","T_CASO4","T_ESP","T_ECE"),"_Count",sep=""))

trops<-aggregate(sum.df["decimal_latitude_corr"], list(sum.df$Taxon),function(x) round(100*(length(x[abs(x)<=23.43780])/length(x)),2))
colnames(trops)<-c("Taxon","TropPcnt")
subtrops<-aggregate(sum.df["decimal_latitude_corr"], list(sum.df$Taxon),function(x) round(100*(length(x[abs(x)>23.43780 & abs(x)<=35.00000])/length(x)),2))
colnames(subtrops)<-c("Taxon","SubTropPcnt")
temp<-aggregate(sum.df["decimal_latitude_corr"], list(sum.df$Taxon),function(x) round(100*(length(x[abs(x)>35.00000 & abs(x)<=60.00000])/length(x)),2))
colnames(temp)<-c("Taxon","TempPcnt")
coord.count<-aggregate(sum.df["decimal_latitude_corr"], list(sum.df$Taxon),function(x) length(x))
colnames(coord.count)<-c("Taxon","NumLatLong")

new.bios<-paste(LETTERS[1:14],"Pcnt",sep="")
for(p in seq(1:14)){
	habs<-aggregate(sum.df["BIOME"], list(sum.df$Taxon),function(x) round(100*(length(x[x==p & !is.na(x)])/sum(!is.na(x))),2))
	colnames(habs)<-c("Taxon",new.bios[p])
	coord.count<-merge(coord.count,habs,by="Taxon")
}

new.plate.codes<-paste(sort(unique(sum.df$Code)),"Pcnt",sep="")
new.plates<-sort(unique(sum.df$Code))
for(p in 1:length(new.plates)){
	plates<-aggregate(sum.df["Code"], list(sum.df$Taxon),function(x) round(100*(length(x[x==new.plates[p] & !is.na(x)])/sum(!is.na(x))),2))
	colnames(plates)<-c("Taxon",new.plate.codes[p])
	coord.count<-merge(coord.count,plates,by="Taxon")
}

sum.table<-merge(meds,counts,by="Taxon")
sum.table<-merge(sum.table,means,by="Taxon")
sum.table<-merge(sum.table,sd.samps,by="Taxon")
sum.table<-merge(sum.table,max.meas,by="Taxon")
sum.table<-merge(sum.table,min.meas,by="Taxon")
sum.table<-merge(sum.table,ranges,by="Taxon")
sum.table<-merge(sum.table,trops,by="Taxon")
sum.table<-merge(sum.table,subtrops,by="Taxon")
sum.table<-merge(sum.table,temp,by="Taxon")
sum.table<-merge(sum.table,coord.count,by="Taxon")

habssies<-c(paste(c("A","B","C","D","E","F","G","H","I","J","K","L","M","N"),"Pcnt",sep=""))
for(x in 1:length(habssies)){
	for(p in 1:nrow(sum.table)){
		if(is.na(sum.table[p,habssies[x]])){
			sum.table[p,habssies[x]]<-0
		}
	}
}
open<-c("GPcnt","HPcnt","IPcnt","JPcnt","KPcnt","MPcnt")
closed<-c("APcnt","BPcnt","CPcnt","DPcnt","EPcnt","NPcnt","FPcnt","LPcnt")

sum.table$Open<-0
sum.table$Closed<-0

for(x in 1:nrow(sum.table)){
	if(sum(sum.table[x,open])>0.33){
		sum.table$Open[x]<-1
	}
	if(sum(sum.table[x,closed])>0.33){
		sum.table$Closed[x]<-1
	}
}

sum.table$habcod<-NA

for(x in 1:nrow(sum.table)){
	if(sum.table$Open[x]==0 & sum.table$Closed[x]==1){
		sum.table$habcod[x]<-0
	}
	if(sum.table$Open[x]==1 & sum.table$Closed[x]==0){
		sum.table$habcod[x]<-1
	}
	if(sum.table$Open[x]==1 & sum.table$Closed[x]==1){
		sum.table$habcod[x]<-"0&1"
	}
}

write.csv(sum.table,file="~/Desktop/ants_ecogeography/ant_summaries_30may2018_wsoils_CLEANED.csv",row.names=FALSE)



#Subsequently noticed 0.33 was used instead of 33.3, so re-did the open/closed delimitation
sum.table<-read.csv(file="~/Documents/papers_reviews/papers/ants_ecogeography/ant_summaries_30may2018_wsoils_CLEANED.csv",stringsAsFactors=FALSE)
open<-c("GPcnt","HPcnt","IPcnt","JPcnt","KPcnt","MPcnt")
closed<-c("APcnt","BPcnt","CPcnt","DPcnt","EPcnt","NPcnt","FPcnt","LPcnt")

sum.table$Open<-0
sum.table$Closed<-0

for(x in 1:nrow(sum.table)){
	if(sum(sum.table[x,open])>33.3){
		sum.table$Open[x]<-1
	}
	if(sum(sum.table[x,closed])>33.3){
		sum.table$Closed[x]<-1
	}
}

sum.table$habcod<-NA

for(x in 1:nrow(sum.table)){
	if(sum.table$Open[x]==0 & sum.table$Closed[x]==1){
		sum.table$habcod[x]<-0
	}
	if(sum.table$Open[x]==1 & sum.table$Closed[x]==0){
		sum.table$habcod[x]<-1
	}
	if(sum.table$Open[x]==1 & sum.table$Closed[x]==1){
		sum.table$habcod[x]<-"0&1"
	}
}

write.csv(sum.table,file="~/Documents/papers_reviews/papers/ants_ecogeography/ant_summaries_30may2018_wsoils_CLEANED_33.3.csv",row.names=FALSE)




#re-did this part because had not used corrected latitude in previous run.  This will only affect STRAPP analyses for latitude, which were re-run on 12 Aug 2018
sum.df<-read.csv(file="~/Desktop/ants_ecogeography/unique_coords_per_species_cleaned_29May2018.wsoils.nona.csv",stringsAsFactors=FALSE)
taxs<-unique(sum.df$Taxon)
columns<-c("Taxon","Occs_BioCl","TropPcnt","SubTropPcnt","TempPcnt","NTempPcnt","STempPcnt","PolarPcnt","NPolarPcnt","SPolarPcnt","SubTempPolPcnt","NSubTempPolPcnt","SSubTempPolPcnt","TropSubPcnt","TempPolPcnt","NTempPolPcnt","STempPolPcnt","SingOcc","SLimit","NLimit","LatRange","MostEquat","MostPolar","MdnAbsLat","MeanAbsLat","SDAbsLat","Elev_Min","Elev_Max","Elev_Range","Elev_Mean","Elev_Mdn","Elev_SD","Elev_Var")
sum.df.sum<-data.frame(matrix(ncol=length(columns),nrow=length(taxs)))
colnames(sum.df.sum)<-columns
sum.df.sum$Taxon<-taxs

for(x in 1:nrow(sum.df.sum)){
	df.dedup.bio<-sum.df[sum.df$Taxon %in% sum.df.sum$Taxon[x],]
	sum.df.sum$Occs_BioCl[x]<-nrow(df.dedup.bio)
	sum.df.sum$TropPcnt[x]<-round(100*(length(df.dedup.bio$decimal_latitude_corr[abs(df.dedup.bio$decimal_latitude_corr)<=23.43780])/sum.df.sum$Occs_BioCl[x]),2)
	sum.df.sum$SubTropPcnt[x]<-round(100*(length(df.dedup.bio$decimal_latitude_corr[abs(df.dedup.bio$decimal_latitude_corr)>23.43780 & abs(df.dedup.bio$decimal_latitude_corr)<=35.00000])/sum.df.sum$Occs_BioCl[x]),2)
	sum.df.sum$TempPcnt[x]<-round(100*(length(df.dedup.bio$decimal_latitude_corr[abs(df.dedup.bio$decimal_latitude_corr)>35.00000 & abs(df.dedup.bio$decimal_latitude_corr)<=60.00000])/sum.df.sum$Occs_BioCl[x]),2)
	sum.df.sum$NTempPcnt[x]<-round(100*(length(df.dedup.bio$decimal_latitude_corr[df.dedup.bio$decimal_latitude_corr>35.00000 & df.dedup.bio$decimal_latitude_corr<=60.00000])/sum.df.sum$Occs_BioCl[x]),2)
	sum.df.sum$STempPcnt[x]<-round(100*(length(df.dedup.bio$decimal_latitude_corr[df.dedup.bio$decimal_latitude_corr<(-35.00000) & df.dedup.bio$decimal_latitude_corr>=-60.00000])/sum.df.sum$Occs_BioCl[x]),2)
	sum.df.sum$PolarPcnt[x]<-round(100*(length(df.dedup.bio$decimal_latitude_corr[abs(df.dedup.bio$decimal_latitude_corr)>60.00000])/sum.df.sum$Occs_BioCl[x]),2)				
	sum.df.sum$NPolarPcnt[x]<-round(100*(length(df.dedup.bio$decimal_latitude_corr[df.dedup.bio$decimal_latitude_corr>60.00000])/sum.df.sum$Occs_BioCl[x]),2)				
	sum.df.sum$SPolarPcnt[x]<-round(100*(length(df.dedup.bio$decimal_latitude_corr[df.dedup.bio$decimal_latitude_corr<(-60.00000)])/sum.df.sum$Occs_BioCl[x]),2)				
	sum.df.sum$SubTempPolPcnt[x]<-round(100*(length(df.dedup.bio$decimal_latitude_corr[abs(df.dedup.bio$decimal_latitude_corr)>23.43780])/sum.df.sum$Occs_BioCl[x]),2)
	sum.df.sum$NSubTempPolPcnt[x]<-round(100*(length(df.dedup.bio$decimal_latitude_corr[df.dedup.bio$decimal_latitude_corr>23.43780])/sum.df.sum$Occs_BioCl[x]),2)
	sum.df.sum$SSubTempPolPcnt[x]<-round(100*(length(df.dedup.bio$decimal_latitude_corr[df.dedup.bio$decimal_latitude_corr<(-23.43780)])/sum.df.sum$Occs_BioCl[x]),2)
	sum.df.sum$TropSubPcnt[x]<-round(100*(length(df.dedup.bio$decimal_latitude_corr[abs(df.dedup.bio$decimal_latitude_corr)<=35.00000])/sum.df.sum$Occs_BioCl[x]),2)
	sum.df.sum$TempPolPcnt[x]<-round(100*(length(df.dedup.bio$decimal_latitude_corr[abs(df.dedup.bio$decimal_latitude_corr)>35.00000])/sum.df.sum$Occs_BioCl[x]),2)
	sum.df.sum$NTempPolPcnt[x]<-round(100*(length(df.dedup.bio$decimal_latitude_corr[df.dedup.bio$decimal_latitude_corr>35.00000])/sum.df.sum$Occs_BioCl[x]),2)
	sum.df.sum$STempPolPcnt[x]<-round(100*(length(df.dedup.bio$decimal_latitude_corr[df.dedup.bio$decimal_latitude_corr<(-35.00000)])/sum.df.sum$Occs_BioCl[x]),2)
	if(nrow(df.dedup.bio)==1){
		sum.df.sum$SingOcc[x]<-"Yes"
		sum.df.sum$MeanAbsLat[x]<-abs(df.dedup.bio$decimal_latitude_corr)
		sum.df.sum$MdnAbsLat[x]<-abs(df.dedup.bio$decimal_latitude_corr)
	}
	if(nrow(df.dedup.bio)>1){
		sum.df.sum$SingOcc[x]<-"No"
		sum.df.sum$SLimit[x]<-min(df.dedup.bio$decimal_latitude_corr)
		sum.df.sum$NLimit[x]<-max(df.dedup.bio$decimal_latitude_corr)
		sum.df.sum$LatRange[x]<-max(dist(df.dedup.bio$decimal_latitude_corr))
		sum.df.sum$MostEquat[x]<-df.dedup.bio$decimal_latitude_corr[which.min(abs(df.dedup.bio$decimal_latitude_corr))]
		sum.df.sum$MostPolar[x]<-df.dedup.bio$decimal_latitude_corr[which.max(abs(df.dedup.bio$decimal_latitude_corr))]
		sum.df.sum$MeanAbsLat[x]<-mean(abs(df.dedup.bio$decimal_latitude_corr))
		sum.df.sum$MdnAbsLat[x]<-median(abs(df.dedup.bio$decimal_latitude_corr))
		sum.df.sum$SDAbsLat[x]<-sd(abs(df.dedup.bio$decimal_latitude_corr))
		sum.df.sum$Elev_Min[x]<-min(df.dedup.bio$ELEV)
		sum.df.sum$Elev_Max[x]<-max(df.dedup.bio$ELEV)
		sum.df.sum$Elev_Range[x]<-max(df.dedup.bio$ELEV)-min(df.dedup.bio$ELEV)
		sum.df.sum$Elev_Mean[x]<-mean(df.dedup.bio$ELEV)
		sum.df.sum$Elev_Mdn[x]<-median(df.dedup.bio$ELEV)
		sum.df.sum$Elev_SD[x]<-sd(df.dedup.bio$ELEV)
		sum.df.sum$Elev_Var[x]<-var(df.dedup.bio$ELEV)
	}
}

write.csv(sum.df.sum,file="~/Desktop/ants_ecogeography/ant_latlong_summaries_30may2018.onlywsoilsandenv.CLEANED.csv",row.names=FALSE)




#THis was re-run, and pct values in ant_summaries_30may2018_wsoils_CLEANED.csv should be fine bc based no corr latitudes, but values in ant_summaries_w_latlong_merged_30may2018_onlywsoilsandenv.CLEANED.csv come from ant_latlong_summaries_30may2018.onlywsoilsandenv.CLEANED.csv, so re-ran this.
sum.table<-read.csv(file="~/Desktop/ants_ecogeography/ant_summaries_30may2018_wsoils_CLEANED.csv", stringsAsFactors=FALSE)
rownames(sum.table)<-sum.table$Taxon
sum.table<-sum.table[,!colnames(sum.table) %in% c("TropPcnt","SubTropPcnt","TempPcnt")]
sum.df.sum<-read.csv(file="~/Desktop/ants_ecogeography/ant_latlong_summaries_30may2018.onlywsoilsandenv.CLEANED.csv", stringsAsFactors=FALSE)
rownames(sum.df.sum)<-sum.df.sum$Taxon
sum.df.sum<-sum.df.sum[,!colnames(sum.df.sum) %in% "Taxon"]
sum.table.combo<-merge(sum.table,sum.df.sum,by=0,all=TRUE)
sum.table.combo<-sum.table.combo[,2:ncol(sum.table.combo)]
colnames(sum.table.combo)[1]<-"Taxon"
write.csv(sum.table.combo,file="~/Desktop/ants_ecogeography/ant_summaries_w_latlong_merged_30may2018_onlywsoilsandenv.CLEANED.csv",row.names=FALSE)
#end of 12 Aug 2018 re-run...only thing changed from previous version used was that pct Trop/Subtrop/Temp based on corrected latitudes...all other values the same and analyses fine.





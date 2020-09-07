library(utils)
library(httr)
library(dplyr)
library("jsonlite")
library(splitstackshape)
library(maditr) 

#download the dataset from the ECDC website to a local temporary file
GET("https://opendata.ecdc.europa.eu/covid19/casedistribution/csv", authenticate(":", ":", type="ntlm"), write_disk(tf <- tempfile(fileext = ".csv")))

#read the Dataset sheet into R
ecdc <- read.csv(tf)

#the country Nambia is coded as NA, hence we need to recode it different
#First make the variable a character (in factor I would've need to add another level)

ecdc=ecdc %>% 
  mutate(
    geo.id=as.character(geoId),
    dateRep=as.Date(dateRep,format="%d/%m/%Y"),
    countries=as.character(countriesAndTerritories))%>%
  tidyr::replace_na(list(geo.id="Nam"))%>%
  select(dateRep:year, geo.id,countries,popData2019,cases,deaths)

#load the data to assign continents

json_file <- 'https://datahub.io/JohnSnowLabs/country-and-continent-codes-list/datapackage.json'
json_data <- fromJSON(paste(readLines(json_file), collapse=""))

# get list of all resources:
print(json_data$resources$name)

# print all tabular data(if exists any)
for(i in 1:length(json_data$resources$datahub$type)){
  if(json_data$resources$datahub$type[i]=='derived/csv'){
    path_to_file = json_data$resources$path[i]
    data <- read.csv(url(path_to_file))
    print(data)
  }
} 

#make arrangement to join the databases
data=data %>% 
  mutate(cont_id=as.character(Continent_Code),
         geo.id=as.character(Two_Letter_Country_Code))%>%
  tidyr::replace_na(list(cont_id="NrA", geo.id="Nam"))%>%
  select(cont_id,Continent_Name,geo.id)

#Join both dataset with the two-letter codes of the countries

sa.cs=plyr::join(data, ecdc, by=c("geo.id"))%>%
 na.omit(cd)%>%
 arrange(Continent_Name, countries, dateRep)%>%
 group_by(geo.id)%>%
 mutate(
  cumsu.c=cumsum(cases),
  cumsu.d=cumsum(deaths))%>%
 filter(cumsu.c!=0)%>%
 getanID(id.vars=c("countries"))%>%#Assing an id that would be a proxy for number of days since the outbreak
 filter(cont_id=="SA")%>%
 filter(countries!="Suriname",countries!="Guyana")%>%
 filter(.id<91)%>%
 filter(geo.id!="FK")%>%
 filter(cases>=0)%>%
 data.frame()
# Eliminate databases to save memory space
rm(ecdc,data) 

####################### TABLES ############################
tab1 = sa.cs%>%
 filter(.id==90)%>%
 mutate(
  inc=round(cumsu.c/popData2019*100000,2),
  mor=round(cumsu.d/popData2019*100000,2))%>%
 select(countries, popData2019, cumsu.c, inc,cumsu.d, mor)%>%
 rename(
  Countries = countries,
  'Total population (day 90)' = popData2019,
  'Cumulative cases (day 90)' = cumsu.c,
  'Incidence rate per 100,000 (day 90)' = inc,
  'Cumulative deaths (day 90)' = cumsu.d,
  'Mortality rate per 100,000 (day 90)' = mor)
write.csv(tab1, "table_1_epi_outcomes.csv")

############ R0 estimation functions ##############
r0.estimate = function(df, days){
 gen <- R0::generation.time(type=c("gamma"), c(mean=3.96,4.75))
 res <- R0::estimate.R(df$cases,GT=gen,methods=c("ML"), end=days)
 return(res)
}

## Base de datos para R0 en diferentes periodos de incubacion, ya que Brazil tuvo
## subreporte en los 1ros dias de la pandemia
sa.cs.b = sa.cs %>% filter(!(countries == "Brazil"&.id<5))
# Base de datos para los cumulative cases
sa.14 = sa.cs%>% filter(.id==14)%>% dplyr::select(countries, cumsu.c)
r0.14 = plyr::dlply(sa.cs,"countries",function(x){r0.estimate(x,13)})
r0.7 = plyr::dlply(sa.cs.b,"countries",function(x){r0.estimate(x,7)})
r0.5 = plyr::dlply(sa.cs.b,"countries",function(x){r0.estimate(x,5)})

r0.composite = function(est){
  r0.co = plyr::ldply(est,function(x){x$estimates$ML$R})%>%rename(R0 = V1)
  r0.cr = plyr::ldply(r0.14,function(x){x$estimates$ML$conf.int})%>%
    rename(lower = V1, upper = V2)
  r0=left_join(r0.co,r0.cr,by=c("countries"))%>%
    mutate(
      R0 = round(R0,2),
      lower = round(lower,2),
      upper = round(upper,2))%>%
    tidyr::unite(col=ci95,lower,upper,sep="-")%>%
    rename('95% CI' = ci95)
  return(r0)}

r0.5.est = r0.composite(r0.5)%>%rename(R0_5d = R0)
r0.7.est = r0.composite(r0.7)%>%rename(R0_7d = R0)
r0 = r0.composite(r0.14)%>%
  left_join(r0.7.est, by = c("countries"), suffix=c("_14d","_7d"))%>%
  left_join(r0.5.est, by = c("countries"))%>%
  left_join(sa.14, by=c("countries"))%>%
  rename(Countries = countries)

############ Rt estimation functions ##############
rt.estimate <- function (df){
 res <- EpiEstim::estimate_R (
   incid = df$cases, 
   method = "parametric_si",
   config = EpiEstim::make_config(list(mean_si = 3.96, std_si = 4.75,
                        t_start = 2:(nrow(df)-5), 
                        t_end = 7:nrow(df))))
 tab <- data.frame(res$R$t_start, res$R$t_end, res$R$Mean, res$R$Quantile.0.025, 
                   res$R$Quantile.0.975)
 names(tab)<-c("day.start","day.end","Rt","lower","upper")
 return(tab)
 }

rt.c=plyr::ddply(sa.cs,"countries", rt.estimate)

rt.90 = plyr::ddply(sa.cs,"countries",function(x){
 rt = rt.estimate(x)
 lt = tail(rt,1)
 return(lt)})

sa.90 = sa.cs%>% filter(.id==90)%>% dplyr::select(countries, cumsu.c)
rt.90 = plyr::join(rt.90,sa.90,by=c("countries"))%>%
 mutate(
  Rt = round(Rt,2),
  lower = round(lower,2),
  upper = round(upper,2))%>%
 tidyr::unite(col=cr95, lower, upper, sep = "-")%>%
 dplyr::select(countries, cumsu.c, Rt, cr95)%>%
 rename(
  Countries = countries,
  'Cumulative cases' = cumsu.c,
  'Rt' = Rt,
  '95% CrI' = cr95)

r0.rt = left_join(r0,rt.90, by=c("Countries")) %>%
  mutate(diff = R0-Rt)
write.csv(r0.rt, "r0_rt90_countries.csv")

####################### FIGURE LOG-SCALE CASES&DEATHS ############################
library(ggplot2)

fig1=sa.cs%>%
 select(countries, .id, cumsu.c, cumsu.d)%>%
 tidyr::gather(key = "status", value = "val", -countries,-.id)%>%
 ggplot(aes(x=.id,y=log10(val),group=status))+
  geom_line(aes(linetype=status), size = .6)+
  ylim("10","100","1k","10k","100k")+
  ylab(expression(paste('Log'[10],' Counts')))+
  scale_x_continuous(breaks=c(25,50,75))+
  xlab("\nDays after first case")+  
  facet_wrap ( ~ countries, ncol = 4)+
  scale_linetype(
   name="Legend",
   labels=c("Cases","Deaths"))+
  theme_bw()+
  theme(legend.position = c(0.75, 0.15))

setwd("~doc")
jpeg("la_counts.jpg", res=300, width = 7, height = 5, units = 'in')
fig1
dev.off ()

################ FIGURE Rt ################
fig2 = ggplot(rt.c, aes(y=Rt,x=day.end))+
  geom_line(size = .4)+
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .25)+
  geom_hline(yintercept = 1, linetype = "dashed", size = .2)+
  facet_wrap(~countries, ncol = 4)+
  xlab("\nDays since first case") +
  labs(y=expression(R[t][]))+
  theme_bw()

jpeg("rt_la.jpg", res=300, width = 7, height = 5, units = 'in')
fig2
dev.off ()
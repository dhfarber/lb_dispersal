#tidyverse version of late blight script
require(ggplot2)
require(dplyr)
require(purrr)
#desktop WD: 'C:/Users/daniel.farber/Documents/WSU/late blight/gradient arrays'
#read data frames
files <-Sys.glob("*.csv")
file_list<-sapply(files,read.csv,header=F,simplify = F, USE.NAMES = T)
#clean the files
#inspect histogram of column sums; fig1c_lr_cleaned is ready! mon6_dg is ready! fig1main is ready! nork ready! 
#p6 ready! shep ready (although doesn't look like a good DG) unknown 26 ready!
uk26_colsums<-data.frame(column = seq(1,96),colsum = colSums(file_list$unknown26.csv))
ggplot(uk26_colsums, aes(x = column,y=colsum)) + geom_point()

wrangle<-function(data) #get data frame of mean incidence by row downwind from focus
  {focus_width<-which.max(apply(data,1,mean)) #returns row of maximum incidence %
  print(paste("Center of focus :",focus_width))
  n<-seq(0,length(data[,1])-focus_width,by=1) #sequence of NS gradient
  mean_inc<-apply(data,1,mean)[focus_width:length(apply(data,1,mean))]
  return(data.frame(n,mean_inc))}
files_wrangled<-sapply(file_list,wrangle,simplify = F, USE.NAMES = T) #create list of all DG dataframes
names(files_wrangled) <- gsub('.csv','',names(files_wrangled)) #remove .csv from titles
plot_names<-c('Field 1 main focus','Field 1 lower right focus','Field 2','Field 3','Field 4','Field 5','Field 6')
field_name_key<-data.frame(Field = plot_names,File = names(files_wrangled))
#write.table(field_name_key,'field_key.txt')
plot_untransformed<-function(df,title)
  {print(ggplot(df,aes(x = n, y = mean_inc)) + geom_point()+labs(title = title,x='Distance from source',y = 'Mean incidence'))}
for(i in 1:length(files_wrangled))
  {plot_untransformed(files_wrangled[[i]],plot_names[i])}
plot_untransformed(files_wrangled[[1]],"trial")
addminval<-function(df)
  {min<-min(df$mean_inc[df$mean_inc>0])/2 #add half of the minimum value to all observations
  df$mean_inc<-df$mean_inc+min
  return(df)}

dg_no0s<-lapply(files_wrangled,addminval)
dg_df<-dg_no0s %>% reduce(left_join,by = 'n')
colnames(dg_df)<-c('n',names(dg_no0s))
#write.csv2(dg_df,'gradients_no0s.csv')
#get a systemically decided-upon c value
#find furthest row n at which mean incidence >=.8
c_val<-c()
for(i in seq(2,length(dg_df)))
  {nas_rm<-dg_df[complete.cases(dg_df[,i]),]
  print(tail(nas_rm$n[nas_rm[,i]>.75],1))
  c_val<-append(c_val,tail(nas_rm$n[nas_rm[,i]>.8],1))} #deal with NAs

#equation format is bad. see https://stackoverflow.com/questions/7549694/adding-regression-line-equation-and-r2-on-graph to fix it
mipd<-function(field,c_val) #fit data to modified inverse power distribution with fixed c_val (should be .5*focus width)
  {l.tot<-log10(field)
  print(paste("length y: ",length(l.tot)))
  l.dist<-log10(field_df$n+c_val)
  print(paste("length x: ",length(l.dist)))
  tmp.mean.bestc<-lm(l.tot~l.dist)
  eq.mean<-paste("y = "," ", round(summary(tmp.mean.bestc)$coefficients[1],2),"*(x+",c_val,")^",
                 round(summary(tmp.mean.bestc)$coefficients[2],2))
  rsq.leg<-paste("R-squared:",round(summary(tmp.mean.bestc)$r.squared,2))
  a<-ggplot(field_df,aes(x = l.dist, y = l.tot)) + geom_point() +
    geom_abline(intercept = summary(tmp.mean.bestc)$coefficients[1],slope = summary(tmp.mean.bestc)$coefficients[2]) + 
    geom_text(x = 2, y = -0.5, label = eq.mean, parse = TRUE)
  #plot(l.tot~l.distbym,ylab=expression("Log"[10]*"(mean lesions per plant)"),
  #     xlab=expression("Log"[10]*"(distance)"),cex.lab = 1.1)
  #abline(summary(tmp.mean.bestc)$coefficients[1],summary(tmp.mean.bestc)$coefficients[2])
  return(list(a,tmp.mean.bestc))}
#MIPD for all fields using c_value where mean incidence =.8 as found above
mipd_list<-c()
for(i in seq(1,7))
  {field_df<-dg_df[complete.cases(dg_df[i+1]),]#data frame no NAs
  mipd_list<-append(mipd_list,mipd(field_df[,i+1],25))}


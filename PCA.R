remove(list=ls())

library(readxl)
library(dplyr)
library(reshape2)
library(ggplot2)
library(shiny)
library(Cairo) 
library(shiny)
library(ggiraph)
library(plotly)
library(htmlwidgets)
library(factoextra)

Chemical_Annotation <- read_excel("Raw data/metabolon_urine.preprocessed.xlsx", 
                                  sheet = "rowData")

Sample_Meta_Data <- read_excel("Raw data/metabolon_urine.preprocessed.xlsx", 
                               sheet = "colData")

assay <- read_excel("Raw data/metabolon_urine.preprocessed.xlsx", 
                    sheet = "assay")



Sample_Data = Sample_Meta_Data
variable_name = "CHEM_ID"

p_data = assay

id_vars = "PARENT_SAMPLE_NAME"
time_points_var  = "TIME_POINT"
Chemical_name_data = Chemical_Annotation

chem_id_variable = "CHEM_ID"
chem_name_variable = "CHEMICAL_NAME"


log_2_fold_change=T


Sample_Data  = Sample_Data %>% group_by( Subject) %>% mutate(order_sample = order(Day))
Sample_Data = Sample_Data %>% mutate(TIME_POINT = paste0(Subject,"_",order_sample))

Sample_Data_2 = ungroup(Sample_Data) %>% select(Sampl,Subject,TIME_POINT,Dose)

Sample_Data_3 = Sample_Data_2 %>% left_join(p_data)

melt_data = melt(Sample_Data_3,id.vars = c(id_vars,time_points_var,"GROUP") , variable.name = variable_name,factorsAsStrings=F)


Chemical_Annotation_2 = Chemical_name_data[,c(chem_id_variable ,chem_name_variable) ]
Chemical_Annotation_2 = Chemical_Annotation_2 %>%
  mutate(across(everything(), as.character))

melt_data_2 = left_join(melt_data,Chemical_Annotation_2)


#melt_data_2$Chemical_time = paste0(melt_data_2[,"CHEMICAL_NAME"],"_", melt_data_2[,time_points_var]) 





melt_data_2[,c(id_vars,variable_name,"GROUP")] = NULL

melt_data_2 = melt_data_2 %>% filter(!CHEMICAL_NAME %in% c("metronidazole", "X-22802"))
#melt_data_2$TIME_POINT = as.numeric(factor(melt_data_2$TIME_POINT))


melt_data_2_scale = melt_data_2 %>% group_by(CHEMICAL_NAME) %>% mutate(log_val = log(value),
                                                                       mean_val = mean(log_val),
                                                                       SD_val = sd(log_val),
                                                                       value = (log_val-mean_val)/SD_val,
                                                                       log_val = NULL,
                                                                       mean_val = NULL,
                                                                       SD_val = NULL
                                                                       )
#cast_data = dcast(melt_data_2,formula = TIME_POINT~CHEMICAL_NAME , mean  )

cast_data = dcast(melt_data_2_scale,formula = TIME_POINT~CHEMICAL_NAME , mean  )



df = cast_data
df[,"TIME_POINT"]=NULL

df_colmean = colMeans(df)

df_mean = matrix(1,nrow=nrow(df),ncol=1)
df_mean = df_mean %*%df_colmean



df2 = df - df_mean
df3= as.matrix(df2)

cov_mat = t(df3) %*% df3



eig_cov_mat = eigen(cov_mat)



eig_vecs = eig_cov_mat$vectors
PCs = df3 %*%  eig_vecs


eig_vals = eig_cov_mat$values

cum_variances = sapply(1:length(eig_vals),function(n){sum(eig_vals[1:n])/sum(eig_vals)  }  )
ind_variances=0
ind_variances[2:length(eig_vals)] =0
ind_variances[1] = cum_variances[1]
ind_variances[2:length(eig_vals)] = sapply(2:length(eig_vals),function(n){cum_variances[n]-cum_variances[n-1]  }  )


loadings_df = as.data.frame(eig_cov_mat$vectors)
loadings_df$CHEMICAL_NAME = colnames(df)
loadings_df_subset = loadings_df %>% select(CHEMICAL_NAME,V1,V2,V3)


graphs_data = melt_data_2_scale %>% left_join(loadings_df_subset)

graphs_data_2 = graphs_data %>% mutate(PC1_point = value*V1,
                                       PC2_point = value*V2,
                                       PC3_point = value*V3
) 

graphs_data_3 = graphs_data_2 %>% group_by(TIME_POINT) %>% summarise(PC1 = sum(PC1_point),
                                                                     PC2 = sum(PC2_point),
                                                                     PC3 = sum(PC3_point))

library(rgl)


time = sapply(graphs_data_3$TIME_POINT,function(x){
  
  strsplit(x,"_")[[1]][1]
})
  
graphs_data_3$time = unname(time)
mycolors <- c('red', 'blue', 'green')
graphs_data_3$color <- mycolors[ as.numeric(as.factor(graphs_data_3$time)) ]

plot3d( 
  x=graphs_data_3$PC1, y=graphs_data_3$PC2, z=graphs_data_3$PC3, 
  col = graphs_data_3$color, 
  type = 's', 
  radius = 2,
  xlab="PC1", ylab="PC2", zlab="PC3")

library(plotly)


fig <- plot_ly(graphs_data_3, x = ~PC1, y = ~PC2, z = ~PC3, color = ~time, colors = mycolors)
fig <- fig %>% add_markers()




htmlwidgets::saveWidget(fig , "3d_PCA_all_SERUM.html")



fig <- plot_ly(graphs_data_3, x = ~PC1, y = ~PC2, color = ~time, colors = mycolors)
fig <- fig %>% add_markers()

htmlwidgets::saveWidget(fig , "2d_PCA_all_SERUM.html")


bar_check = loadings_df_subset %>% mutate( abs_pc1= abs(V1),
                                        rank_pc1 = rank(-abs_pc1) )

loadings_bar_check = bar_check %>% filter(rank_pc1<=20)


loadings_p <- plot_ly(loadings_bar_check,x = ~V1, y = ~reorder(CHEMICAL_NAME, V1), type = 'bar', orientation = 'h')





loadings_p = loadings_p %>% layout(
       xaxis = list(title = "Loading"),
       yaxis = list(title =""))

#loadings_p%>% layout(yaxis=list(tickvals=~CHEMICAL_NAME, ticktext=~CHEMICAL_NAME))


htmlwidgets::saveWidget(loadings_p , "PC1_all_loading_SERUM.html")

ind_df = tibble(Variance = round(ind_variances*100,1), component = paste("PC",1:length(ind_variances)))

ind_df = tibble(Variance = round(ind_variances*100,1), component = paste("PC",1:length(ind_variances)))


scree <- plot_ly(ind_df[1:5,], x = ~component, y = ~Variance, type = 'scatter', mode = 'lines+markers') %>% 
  layout(yaxis = list(title = "Variance (%)"),xaxis = list(title = "Component"))

htmlwidgets::saveWidget(scree , "scree_all_SERUM.html")



write.csv(graphs_data_3, "C:/Users/Rik/OneDrive/Metobolites/PCA/PCA_coords_SERUM_all.csv",row.names = F)


write.csv(bar_check, "C:/Users/Rik/OneDrive/Metobolites/PCA/loadings_SERUM_all.csv",row.names = F)

write.csv(ind_df, "C:/Users/Rik/OneDrive/Metobolites/PCA/variances_SERUM_all.csv",row.names = F)
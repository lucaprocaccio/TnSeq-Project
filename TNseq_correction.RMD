---
title: "TNseq_correction"
output: 
  rmdformats::robobook:
    self_contained: TRUE
    highlight: kate
---

```{r setup, include=FALSE}
library(knitr)
library(rmdformats)

## Global options
## Global options
options(max.print="75")
opts_chunk$set(cache=TRUE,
               prompt=FALSE,
               tidy=TRUE,
               comment=NA,
               message=FALSE,
               warning=FALSE)
opts_knit$set(width=75)
```

Name: Luca

Surname: Procaccio

ID_number: 965948

# Introduction

Transposon sequencing (Tn-seq) is a sound method for studying genes and their interactions in microorganisms. This approach is based on the assembly of a saturated transposon insertion library. The library grews under a test condition: changes in frequency of each insertion mutant are determined by sequencing the flanking regions en masse. It follows that these changes reflect the effect of insertion on fitness and they are a quantitative measure of growth rate.
Singolar fitness values for each gene are obtained by averaging over all insertions in that specific gene.

According to the results, genes can be divided into four categories:

- neutral (fitness = 0.96–1.04)
- advantageous (fitness > 1.04)
- disadvantageous (fitness < 0.96)
- possibly essential (fitness = 0)

Genes in the latter category correspond to those in which transposon insertions were absent in the sequenced library. Mutations in these genes either prevented or substantially slowed bacterial replication during outgrowth of the transposon insertion strains.
The other categories take into account detrimental insertions, which decrease in frequency lowering the fitness, advantageous insertions, which increase in frequency increasing the fitness. In neutral category, the number of insertions do not alter the fitness.

# Data

Data come from the work of van Opijnen, Bodi, Camilli (2009), *Tn-seq: high-throughput parallel sequencing for fitness and genetic interaction studies in microorganisms* (https://doi.org/10.1038/nmeth.1377).
In this study the Tn-seq approach was used on a S. pneumoniae culture. Specifically, this culture was transformed with the magellan6 transposon which contains an MmeI restriction site. The result is a bacterial pool in which each bacterium contains a single transposon insertion in its genome. The library is grown in standard broth medium: during the growth detrimental insertions deacreased in frequency whereas advantageous ones increased. Then fitness was computed as written in the introduction.

First, the library useful for the project is attached.

```{r warning=FALSE, message=FALSE}
library(data.table)
```

Data are loaded.

```{r}
# read the txt files and transform them into data.frame
fData<-read.delim("Tn_seq_fitness_data_Opijnen_et_al_2009.txt",header=TRUE,stringsAsFactors = FALSE,sep = "\t")
geneCoord<-read.delim("GCF_000006885.1_ASM688v1_genomic_olt.txt",header=FALSE,stringsAsFactors = FALSE,sep = "\t")

# rename the columns of geneCoord
geneCoord <- reshape::rename(geneCoord, c(V1 = "Locus", V2 = "Start", V3 = "End"))
```

Structure of the two dataframes.

```{r}
head(fData)
```

```{r}
head(geneCoord)
```

Data manipulation, in order to have the same loci in the two dataframes and remove NA values.

```{r}
# create a vector of the loci both present in fData and geneCoord
m<-match(fData$locus,geneCoord$Locus)

# create a vector with NA values present in vector m
w<-which(is.na(m))

# remove NA values from fData
fData<- fData[-w,]

# remove NA values from m
m<-m[-w]

# get only coordinates without NA values in geneCoord
geneCoord<-geneCoord[m,]
```

# Purpose of the project

Since DNA replication induces a bias in fitness values, the purpose of this project is to correct that error. In fact, near the replication fork, the amount of DNA found is present in a larger quantity than that of other areas in the genome: these unbalanced quantities cause an overestimation of the fitness near replication fork. Knowing that the bacterial genome is circular, these regions are at the edges of the linearized chromosome. The results are a "smile" shape when the fitness of a gene is plotted against its position in the genome and a misclassification of some genes.

The graph below shows the "smile" trend for all the data.

```{r}
plot((geneCoord$Start[fData$average_fitness!=0]+geneCoord$End[fData$average_fitness!=0])/2, fData$average_fitness[fData$average_fitness!=0], ylab = "gene_fitness",xlab="gene_middle_coordinate", main="fitness vs genome location")
line_parameters_0 <- loess.smooth((geneCoord$Start[fData$average_fitness!=0]+geneCoord$End[fData$average_fitness!=0])/2, fData$average_fitness[fData$average_fitness!=0])
lines(line_parameters_0, col = "red", lwd = 2)
```

# Methods of work 

In the assignment were suggested two type of approaches:

- Working on windows of 100000 nucleotides and computing the average fitness of the these windows. 
- Working on group of 100-200 genes according to their order and computing the average fitness of these groups.  

The idea is that the vast majority of the genes should be in the neutral category, so it is predictable that all the windows and groups will have similar average fitness values.

I have decided to use both the approaches.

## First method: window

```{r}
# set the size length of the window
size <- 100000

# create a unique dataframe with all the data needed, excluding data with fitness equal to zero
df <- data.frame(locus = fData$locus[fData$average_fitness!=0],
                 start = geneCoord$Start[fData$average_fitness!=0],
                 end = geneCoord$End[fData$average_fitness!=0],
                 fitness = fData$average_fitness[fData$average_fitness!=0],
                 interval = cut(geneCoord$Start[fData$average_fitness!=0], 
                 breaks = seq(0,max(geneCoord$End[fData$average_fitness!=0]) + size,    size)))

# convert data frame to data table 
setDT(df)

# find average fitness for each interval 
interval_average_fitness <- df[ ,list(mean=mean(fitness)), by=interval]
# take the maximum coordinate for each interval 
a <- df[ ,list(end = max(end)), by= interval]
# take the minimum coordinate for each interval
b <- df[ ,list(start = min(start)), by= interval]
# count the number of times each interval is repeated in the dataframe
count <- as.vector(table(df$interval))

# compute middle coordinate for each interval and approximate
middle_coordinate <- round((a$end - b$start)/2, digits = 0)
# find the real middle coordinate for each interval
real_middle_coordinate <- b$start + middle_coordinate

# compute the middle coordinate for each gene
df$gene_middle_coordinate <- round((df$end + df$start)/2, 0)

# update the dataframe with the new values found
df$window_middle_coordinate <- rep(real_middle_coordinate, count)
df$window_fitness <- rep(interval_average_fitness$mean, count)

# add labels according to the average fitness of intervals 
df$window_label <- cut(df$window_fitness, breaks = c(0, 0.96, 1.04, Inf), 
                       labels = c("disadvantageous","neutral","advantageous"))

# add labels according to the fitness of genes 
df$gene_label <- cut(df$fitness, breaks = c(0, 0.96, 1.04, Inf), 
                     labels = c("disadvantageous","neutral","advantageous"))
```

Structure of the resulting dataframe
```{r}
head(df)
```

The average fitness for each interval is shown in the graph below.
As expected all the windows are in the neutral category. Even if all the windows are in the same category the "smile" shape is clearly visible also here.

```{r}
plot(df$window_middle_coordinate, df$window_fitness, col = "red", lwd = 2, xlab = "window_middle_coordinate", ylab = "window_average_fitness")
legend("bottomright", "neutral", col = "black", cex = 0.7, fill = "red")
```


## Second method: group of genes

```{r}
# set the group names
group_names <- c("group_1","group_2","group_3","group_4","group_5",
                 "group_6","group_7","group_8")

# create a vector with the values for which each group name has to be repeated.
repetition <- c(rep(200,nrow(df)%/%200),nrow(df)%%200)

# create a new unique dataframe
new_df <- data.frame(locus = df$locus,
                     start = df$start,
                     end = df$end,
                     fitness = df$fitness,
                     gene_group = rep(group_names,repetition),
                     gene_label = df$gene_label)

#convert dataframe to data table 
setDT(new_df)

# find average fitness for each group 
group_average_fitness <- new_df[ ,list(mean=mean(fitness)), by=gene_group]
# take the maximum coordinate for each group
c <- new_df[ ,list(end = max(end)), by= gene_group]
# take the minimum coordinate for each group
d <- new_df[ ,list(start = min(start)), by= gene_group]

# compute middle coordinate for each group and approximate
middle_coordinate_2 <- round((c$end - d$start)/2, digits = 0)
# # find the real middle coordinate for each group
real_middle_coordinate_2 <- d$start + middle_coordinate_2

# update also this new dataframe with the gene middle coordinate
new_df$gene_middle_coordinate <- df$gene_middle_coordinate

# update the new dataframe with the new values found
new_df$group_middle_coordinate <- rep(real_middle_coordinate_2, repetition)
new_df$group_fitness <- rep(group_average_fitness$mean, repetition)

# add labels according to the average fitness of groups 
new_df$group_label <- cut(new_df$group_fitness, breaks = c(0, 0.96, 1.04, Inf), 
                       labels = c("disadvantageous","neutral","advantageous"))

```

Structure of the resulting new dataframe.

```{r}
head(new_df)
```

The average fitness for each group is shown in the graph below.
Also here all the groups are in the neutral category and the "smile" shape is still present.

```{r}
plot(new_df$group_middle_coordinate, new_df$group_fitness, col = "red", lwd = 2, xlab = "group_middle_coordinate", ylab = "group_average_fitness")
legend("bottomright", "neutral", col = "black", cex = 0.7, fill = "red")
```

Comparison between the elements found in each category for original data, window method and group method.

```{r}
table(df$gene_label)
```

```{r}
table(df$window_label)
```

```{r}
table(new_df$group_label)
```

As it is already written, the totality of windows and groups average fitness values fall in neutral category, while for original data are present elements also in the other two categories.

```{r}
par(mfrow = c(1,3))
barplot(table(df$window_label), col = c("black","red", "green"), main = "Window method", cex.names = 0.9)
barplot(table(df$gene_label), col = c("black","red", "green"), main = "Original", cex.names = 0.9)
barplot(table(new_df$group_label), col = c("black","red", "green"), main = "Group method", cex.names = 0.9)
invisible(dev.off())
```

Graph of original data that accounts also for the category of genes.

```{r}
plot((df$end - df$start)/2 +df$start, df$fitness, col = df$gene_label, xlab = "gene_middle_coordinate", ylab = "gene_fitness")
legend("bottomright",names(table(df$gene_label)), col = c("black","red","green"), cex =0.7, fill = c("black","red","green"))
line_parameters <- loess.smooth((df$end - df$start)/2 +df$start, df$fitness)
lines(line_parameters, col = "blue", lwd = 2)
```

# Correction

To correct the bias introduced by the active replication of the chromosome, it must be used a linear regression model. But, since the genome is circular and replicated in both the directions there is the necessity to transform linear coordinates into circular ones. This can be done by mean of radians.

## Window model

After the transformation of linear coordinates into radians, data are split randomly into train and test set.

```{r}
# transform data into radians
window_radians <- (df$window_middle_coordinate *2*pi) / max(geneCoord$End)
gene_radians <- df$gene_middle_coordinate * 2*pi / max(geneCoord$End)

# create a copy of the first dataframe
df_with_radians <- df
# add coordinates in radians to the copy
df_with_radians$window_radians <- window_radians
df_with_radians$gene_radians <- gene_radians

# create train and test set
train_indices <- sort(sample(1:nrow(df_with_radians), (nrow(df_with_radians)/2) + 1))
train_window <- df_with_radians[train_indices,]
test_window <- df_with_radians[-train_indices,]
```

It is now possible to create a linear regression model. To improve the fitting it is used cosine function within the equation.

```{r}
# create the first model 
first_model <-lm(window_fitness ~ cos(window_radians), data=train_window)
summary(first_model)
```

What it is important to say is that by looking at the adjusted R-squared it can be assumed that the model fits nicely. Both the coefficients have a low p-value, so they are statistically significant. Finally, the p-value associated to the F-statistic is really low, this means that the regression model fits the data better than a model with no predictors.

Homoscedasticity is a topic to keep into account in a linear regression model. It describes a situation in which the error term is the same across all values of the independent variables. The plot below seems to indicate that the residuals and the fitted values are uncorrelated, so the model it is not biased.

```{r}
plot(first_model$fitted.values, first_model$residuals, xlab = "window_model_fitted.values", ylab = "window_model_residuals")
```

Once the model is trained, it must be used on the test set, to verify whether it fits well on new data or if it is present overfitting in the training.
As a measure of goodness of fit it is used the mean squared error(mse). In a model with no errors the mse is equal to zero.
In this case the mse is really low, so the model fits well on test data.

```{r}
# fit the model on test set
prediction_window <- predict(first_model, test_window, se.fit = TRUE)
# compute the mean squared error
mse_1 <- 1/nrow(test_window)*sum(prediction_window$fit - test_window$window_fitness)^2
mse_1
```

Also graphically it can be seen that the model fits quite well.

```{r}
plot(test_window$window_middle_coordinate,test_window$window_fitness, ylab = "test_window_fitness", xlab = "test_window_middle_coordinate")
lines(test_window$window_middle_coordinate,prediction_window$fit, col="red")
```

Now the model is fitted on original data to correct the "smile" shape and to adjust for incorrect classification of genes. First of all, the coefficents of the regression model are obtained and used in the liner model equation using, in this case, as x the middle coordinates in radians for every single gene. In this way are estimated the predicted fitness values. Then the corrected fitness values are computed by subtracting the original fitness values to the predicted ones and summed to one, since they are distributed around this value. Finally, the labels are re-assigned for the new fitness values.

```{r}
fitness_predicted_1 <- first_model$coefficients[1] + first_model$coefficients[2]*cos(df_with_radians$gene_radians)
fitness_corrected_1 <- df_with_radians$fitness - fitness_predicted_1 + 1
label_corrected_1 <- cut(fitness_corrected_1, breaks = c(0, 0.96, 1.04, Inf), 
                         labels = c("disadvantageous","neutral","advantageous"))
```

The graph below shows how the "smile" trend has been corrected.

```{r}
line_parameters_1 <- loess.smooth(df$gene_middle_coordinate, fitness_corrected_1)
plot(df$gene_middle_coordinate, fitness_corrected_1, col = label_corrected_1, xlab = "gene_middle_coordinate")
lines(line_parameters_1, col = "blue", lwd = 2)
legend("bottomright",names(table(df$gene_label)), col = c("black","red","green"), cex =0.7, fill = c("black","red","green"))
```

Comparison between the elements assigned to the three categories in original data and after the correction with window model.

```{r}
table(label_corrected_1)
```

```{r}
table(df$gene_label)
```

Selection of the genes that have changed the label and analysis of their distribution.
It can be noticed that the majority of the changes in label has been made both in the first and last genes in the linearized genome, corresponding to the genes near replication fork. Specifically there has been a meaningful reduction in advantageous category. This is in accordance with the fact that there is an overestimation of fitness in the areas near the replication fork.

```{r}
changing_genes_1 <- which(df$gene_label != label_corrected_1)
hist(changing_genes_1, freq = FALSE, breaks = 100, col = "orange",
     main = "Distribution of changing genes in window model")
```

## Group model

The steps are the same of window model.

```{r}
# create a copy of the second dataframe and transform data into radians
new_df_with_radians <- new_df
group_radians <- (new_df$group_middle_coordinate *2*pi) / max(geneCoord$End)

# add coordinates in radians to the copy 
new_df_with_radians$group_radians <- group_radians
new_df_with_radians$gene_radians <- gene_radians

# create train and test set
train_indices_2 <- sort(sample(1:nrow(new_df_with_radians), (nrow(new_df_with_radians)/2) + 1))
train_group <- new_df_with_radians[train_indices_2,]
test_group <- new_df_with_radians[-train_indices_2,]
```

It is created a linear regression model using also here the periodic function cosine.

```{r}
# create the second model
second_model <-lm(group_fitness ~ cos(group_radians), data=train_group)
summary(second_model)
```

Also in this case the output of the trained regression model is significant, the adjusted R-squared is even higher than the window model, the coefficients are significant having low p-values and the F-statistic is good.

Homoscedasticity seems to be verified because of the absence of an evident pattern, anyway it's not very clear since 8 groups are not that much.

```{r}
plot(second_model$fitted.values, second_model$residuals)
```

The model is used on the test set and the mean squared error is computed.

```{r}
# fit the model on test set
prediction_group <- predict(second_model, test_group, se.fit = TRUE)
# compute mean squared error
mse_2 <- 1/nrow(test_group)*sum(prediction_group$fit - test_group$group_fitness)^2
mse_2
```

Graphic view of the goodness of fit of the group model on test set.

```{r}
plot(test_group$group_middle_coordinate,test_group$group_fitness, xlab = "test_group_middle_coordinate", ylab = "test_group_fitness")
lines(test_group$group_middle_coordinate,prediction_group$fit, col="red")
```

Fitting of the second model on the totality of data with the same procedure already seen for the window model and adjustment for incorrect label classification.

```{r}
fitness_predicted_2 <- second_model$coefficients[1] + second_model$coefficients[2]*cos(new_df_with_radians$gene_radians)
fitness_corrected_2 <- new_df_with_radians$fitness - fitness_predicted_2 + 1
label_corrected_2 <- cut(fitness_corrected_2, breaks = c(0, 0.96, 1.04, Inf), 
                         labels = c("disadvantageous","neutral","advantageous"))
```

Graph after the correction with group model: the smile has been linearized.

```{r}
line_parameters_2 <- loess.smooth(new_df$gene_middle_coordinate, fitness_corrected_2)
plot(new_df$gene_middle_coordinate, fitness_corrected_2, col = label_corrected_2, xlab = "gene_middle_coordinate")
lines(line_parameters_2, col = "blue", lwd = 2)
legend("bottomright",names(table(df$gene_label)), col = c("black","red","green"), cex = 0.7, fill = c("black","red","green"))
```

Comparison between the elements assigned to the three categories in original data and after the correction with group model.

```{r}
table(label_corrected_2)
```

```{r}
table(new_df$gene_label)
```

Selection of the genes that have changed their label and analysis of their distribution.
The same conclusions for window model are valid also for group model.

```{r}
changing_genes_2 <- which(df$gene_label != label_corrected_1)
hist(changing_genes_2, freq = FALSE, breaks = 100, col = "orange", 
     main = "Distribution of changing genes in group model")
```

It is possible to notice that all the genes that change their label with window model are the same genes that change their label with group model.

```{r}
table(changing_genes_1 == changing_genes_2)
```

# Conclusions

The correction was performed following two different approaches, both valid and well performing. Specifically, the gene group model has a higher adjusted R-squared and a higher mean squared error, the window model has a lower adjusted R-squared and a lower mse. However, the correction is good for both the models and the outputs for the corrected labels are the same. 

















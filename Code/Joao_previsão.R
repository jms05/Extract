col = terrain.colors(length(pData(data)$sample), alpha = 0.6)
dataBox = log2(exp)



boxplot(dataBox, col = col, whiskcol = col, pch = 21, cex = 0.8, xlab = 'Amostras', ylab = 'Intensidades (log2)', notch = TRUE, frame = FALSE)


sum(!complete.cases(log2(exp_norm)))

data_antes = data[,data$time=="baseline"]
data_depois = data[,data$time=="1 month"]


addmargins(table(data_antes$agent,data_antes$other))
addmargins(table(data_depois$agent,data_depois$other))




dim(data)
sum(!complete.cases(exp_filtrado))

length(featureNames(data))

data

arrayString = c("Olam","Mundo")
length(arrayString)
arrayString[1]


data


valoresOk = rownames(exp_filtrado)

 
length(valoresOk)

dim(exp_filtrado2)
class(valoresOk)


getFromExpression <- function(expressionFonte,valuesOK)
{
  j=1
  nomes = featureNames(expressionFonte)
  indicesValuesOK =c()
  for(i in 1:length(valuesOK)){
    nomeOK = valuesOK[i]
    indices = which( nomes== nomeOK)
    if(length(indices) > 0){
      indicesValuesOK[j] = indices

      j=j+1
    }
    
    
  }
  return(expressionFonte[indicesValuesOK,])
}

library('genefilter')
exp_filtrado2=data[sds>=2*m,]
exp_filtrado2_matrix = exprs(exp_filtrado2)

exp_filtrado2 = getFromExpression(exp_filtrado2,rownames(exp_filtrado2_matrix[complete.cases(exp_filtrado2_matrix),]))
biocLite("affydata")
biocLite("affy")

library(affydata)
library(affy)
normalized_exp_filtrado2=affy.scalevalue.exprSet(exp_filtrado2,sc=1)
normalized_exp_filtrado2_matrix = exprs(normalized_exp_filtrado2)
c(mean(normalized_exp_filtrado2_matrix),sd(normalized_exp_filtrado2_matrix))

summary(exp_filtrado2)
complete.cases(exp_filtrado2)

if (require(affydata)) {
  data(Dilution)
  normalize.AffyBatch.scaling(Dilution)
}

install.packages("imager")
library(imager)
im<-load.image("Rplot.png")
plot(im,axes=FALSE)

biocLite("ExpressionNormalizationWorkflow")
library(ExpressionNormalizationWorkflow)
normalized_exp_filtrado2=affy.scalevalue.exprSet(exp_filtrado2,sc=1)
normalized_exp_filtrado2_matrix = (exprs(normalized_exp_filtrado2))
c(sd(normalized_exp_filtrado2_matrix),mean(normalized_exp_filtrado2_matrix))



##Divis~ão do dataset para treino e teste
percentDataTrain = 0.75
nomesComplete =rownames(exp_filtrado2_matrix[complete.cases(exp_filtrado2_matrix),])
nomesComplete_shufle= sample(nomesComplete)
dataToTrain = as.integer(percentDataTrain*length(nomesComplete_shufle))

exp_filtrado2_test =getFromExpression(exp_filtrado2,nomesComplete_shufle[1:dataToTrain])
exp_filtrado2_train = getFromExpression(exp_filtrado2,nomesComplete_shufle[1:dataToTrain])

slotNames(exp_filtrado2_test)
exp_filtrado2_test$other



library(knn)
valores.previstos=knn(t(exprs(exp_filtrado2_test)),exprs(exp_filtrado2_test)[, 1],exp_filtrado2_train$other)
exprs(exp_filtrado2_test)[, 1]
valores.previstos
table(valores.previstos,exp_filtrado2_test$other)
dim(exp_filtrado2_train)

df <- data.frame(matrix(geneNames, nrow=132, byrow=T))
matrix(geneNames)

normalized_exp_filtrado2_matrix


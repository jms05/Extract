########## Relatório ############ 

#******************carregar dados e visualizar os dados existentes**************************************
install.packages("gplots")
install.packages("dendextend")#para por por no heatmap nas linhas
install.packages("dendextendRcpp")#para por por no heatmap nas linhas
install.packages("d3heatmap")#para fazer o heatmap
source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("GEOquery")
biocLite("limma")
biocLite("genefilter")

biocLite("hgu133a.db")
library("hgu133a.db")
#biocLite("hgu95av2.db")

biocLite("illuminaHumanv4.db")

library("illuminaHumanv4.db") #Get this library if you don't have

biocLite("hgu133plus2.db")
library("hgu133plus2.db")

library(Biobase) 
library(GEOquery)
library("hgu95av2.db")
library('genefilter')
library(limma)
library(gplots)
library(dendextend)#heatmap
library(dendextendRcpp)#heatmap
library(d3heatmap)#heatmap




#download do ficheiro GDS
gds <- getGEO('GDS5393', destdir=".")  
gds
#converte GDS em ExpressionSet e ao mesmo tempo faz a tranforma??o logar?tmica
data <- GDS2eSet(gds,do.log2=TRUE)  
dim(data)#Features:48107 Samples:120
exp = exprs(data)#matriz num?rica	(genes nas linhas, amostras	nas	colunas)
exp
colnames(exp)#amostras 
rownames(exp)#nomes genes
class(exp)#matriz
sampleNames(data)#nome das samples
featureNames(data)#nome das features
varMetadata(data)#atributos dos metadados das amostras
#n?o temos nenhuma descricao e os atribitos sao 5: sample;agent, other, time,individual;description
data$sample[1:5]#120 levels que correspondem as amostras
data$agent[1:5]#2 levels:control lithium
data$other[1:5]#2 levels: non-responder responder
data$individual[1:5]#60 levels que correspondem aos pacientes
data$time #1 mes e baseline 

exp
data$description[1:5]# [1] "Value for GSM1105438: OPT_1-baseline + Li+OPT_1 ; src: Peripheral blood" 
annotation(data)# VER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
experimentData(data)
abstract(data)# "Analysis of peripheral blood from patients with bipolar disorder 
#before and 1 month after lithium treatment. Response of patients to lithium assessed 
#after 6 months. Results identify a gene expression
#signature for the response to lithium treatment in patients with bipolar disorder."

#**********************************PrÈ-processamento dos dados**********************************************

#-> tratar valores omissos
sum(is.na(exp))
#Foram encontrados 94080 valores omissos no nosso DataSet. De forma a tornar possÌvel a an·lise dos dados, procedeu-se
#‡ remoÁ„o dos valores omissos. 
exp_complete=na.exclude(exp)

#-> normalizacao
#Para verificar se os nossos dados seguiam uma distribuiÁ„o normal, calcul·mos a mÈdia e desvio padr„o dos dados do nosso 
#DataSet. 
mean(exp_complete) #2.88391
sd(exp_complete)#0.2084394
#Verific·mos que os nossos dados n„o seguiam uma distribuiÁ„o normal, pelo que se procedeu ‡ normalizaÁ„o dos mesmos. 

exp_norm=scale(exp_complete) 
dim(exp_norm)
mean(exp_norm)#5.533452e-17
sd(exp_norm)#0.9999895


#De seguida, procurou-se remover os genes que possuiam valores pouco dispares ao longo das amostras, uma vez que n„o 
#apresentam variabilidade significante para a distinÁ„o dos grupos amostrais. 
#Para tal, procedeu-se ‡ seleÁ„o dos genes cujo desvio padr„o dos valores de express„o ao longo das amostras È
#duas vezes maior do que a mediana dos desvios padr„o de todos os genes.
#Com este passo, reduziu-se os dados para posterior an·lise para 7588 features. 


sds=rowSds(exp_norm)  #desvio padrao de cada linha 
m=median(sds)
hist(sds,breaks=50,col='mistyrose')
abline(v=m*2,col="red",lwd=4,lty=2)#usamos este
exp_filtrado=exp_complete[sds>=2*m,]
exp_filtrado #7588 features, 64 samples (utilizado nas an?lises seguintes)


#para construir o histograma tivemos que usar a matriz 

#Histograma filtrado
sds__1=rowSds(exp_filtrado)
hist(sds__1,breaks=50,col='mistyrose')

####### Analise de expressao diferencial 
library('genefilter')
exp_filtrado2=data[sds>=2*m,]
#para a analise diferencial temos que usar o expression set dai fazer a linha anterior


#######Análise de expressão diferencial e Enriquecimento#######
biocLite("affydata")
biocLite("affy")
<<<<<<< HEAD
biocLite("affyPLM")
library(affydata)
library(affy)
library(affyPLM)
=======
library(affydata)
library(affy)

>>>>>>> 0e15e90cf70421776e1fa0ec23f5e0eb4a6c4f07
###### Para realizarmos testes de hipóteses paramétricos (t-test) precisamos
#de normalizar o nosso ExpressionSet, mantendo os filtros (remoção de NA e aplicação de filtros "Flat Patterns")
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

exp_filtrado2_matrix = exprs(exp_filtrado2)
exp_filtrado2 = getFromExpression(exp_filtrado2,rownames(exp_filtrado2_matrix[complete.cases(exp_filtrado2_matrix),]))
<<<<<<< HEAD
exp_filtrado2_matrix = scale(exprs(exp_filtrado2))
sum(is.na(exp_filtrado2_matrix))
mean(exp_filtrado2_matrix)
sd(exp_filtrado2_matrix)
phenoData1=phenoData(exp_filtrado2)
featData1=featureData(exp_filtrado2)
expData1=experimentData(exp_filtrado2)
normalized_exp_filtrado2=ExpressionSet(assayData = exp_filtrado2_matrix, phenoData = phenoData1, featureData = featData1, experimentData = expData1)
normalized_exp_filtrado2
mat=exprs(normalized_exp_filtrado2)
mean(mat)
sd(mat)

#####Testes estatísticos

#Neste passo foram identificados, através de testes de hipóteses, os 25 genes 
#que apresentavam menor p-value, ou seja, os genes cuja evidência estatística da
#expressão génica diferencial foi mais significativa. 
#Os resultados obtidos foram registados nas seguintes tabelas e organizados segundo
#o nome do gene, identificador do gene e respetiva função. 

#########Controlo x Lítio (1 mês de tratamento) ############

#Na seguinte tabela encontram-se os 25 genes com maior evidência estatística de expressão génica
#diferencial no que diz respeito ao atributo "agent". Este teste foi feito com o intuito
#de perceber se, após 1 mês de tratamento em questão (com lítio), existem diferenças estatisticamente
#significativas na expressão dos genes. 
#De facto, constatou-se que respectivos p-values (p-values < 0.01) permitem rejeitar a hipótese nula, pelo
#que podemos assumir que os genes das amostras em condição de controlo e condição de tratamento
#com lítio são diferencialmente expressos. Numa primeira análise, isto pode indiciar 
#uma influencia do lítio na expressão destes genes.  


test1=rowttests(normalized_exp_filtrado2[normalized_exp_filtrado2$time=="1 month"],"agent")
test1
sign_pvalue1=test1[test1$p.value<0.01,]
pvalues1Signif=sign_pvalue1[order(sign_pvalue1$p.value),]
nomes1=rownames(pvalues1Signif)[1:25]
geneNames1=unlist(mget(nomes1, illuminaHumanv4SYMBOL))
geneNames1
geneFunction1=unlist(mget(nomes1, illuminaHumanv4GENENAME))
geneFunction1


#######Controlo x Litio (baseline)#########

#Na seguinte tabela encontram-se os 25 genes com maior evidência estatística de expressão génica
#diferencial no que diz respeito ao atributo "agent". Este teste foi feito com o intuito
#de perceber se existem diferenças estatisticamente significativas na expressão dos genes no início da experiência. 
#De facto, constatou-se que respectivos p-values (p-values < 0.01) permitem rejeitar a hipótese nula, pelo
#que podemos assumir que os genes das amostras em condição de controlo e condição de tratamento com lítio
#no ínicio do estudo são diferencialmente expressos. Numa primeira análise, isto pode indiciar 
#uma influencia do lítio na expressão destes genes. 
#Tendo em conta os genes que obtivemos na tabela anterior, relativos a 1 mês de tratamento, podemos verificar
#que estes não são os mesmos obtidos na tabela seguinte. Estes resultados suportam a hipótese de que os genes obtidos na
#tabela anterior possam ser, de facto, influenciados ao nível da sua expressão pelo tratamento com lítio, uma vez que estes 
#genes não apresentavam a maior diferença de expressão, comparativamente a outros genes, no início da experiência. 

test2=rowttests(normalized_exp_filtrado2[normalized_exp_filtrado2$time=="baseline"],"agent")
test2
sign_pvalue2=test2[test2$p.value<0.01,]
pvalues2Signif=sign_pvalue2[order(sign_pvalue2$p.value),]
nomes2=rownames(pvalues2Signif)[1:25]
geneNames2=unlist(mget(nomes2, illuminaHumanv4SYMBOL))
geneFunction2=unlist(mget(nomes2, illuminaHumanv4GENENAME))
geneFunction2
#unique(as.data.frame(table(nomes2, geneNames2))[1])
#unique(as.data.frame(table(nomes2, geneNames2))[2])
#unique(as.data.frame(table(nomes2, geneFunction2))[2])
#c(unique(as.data.frame(table(nomes2, geneNames2))[2]), unique(as.data.frame(table(nomes1, geneNames1))[2]))
=======
normalized_exp_filtrado2=affy.scalevalue.exprSet(exp_filtrado2,sc=1)
normalized_exp_filtrado2_matrix = exprs(normalized_exp_filtrado2)
c(mean(normalized_exp_filtrado2_matrix),sd(normalized_exp_filtrado2_matrix))



#########Controlo x Lítio (1 mês de tratamento) ############ 
test1=rowttests(normalized_exp_filtrado2[normalized_exp_filtrado2$time=="1 month"],"agent")
test1
pvalues1= test1$p.value
pvalues1
pvalues1Signif= order(test1[test1$p.value<0.01, "p.value"])
pvalues1Signif
indices_melhoresGenes=pvalues1Signif[1:25]
indices_melhoresGenes
geneFeatNames=featureNames(normalized_exp_filtrado2[indices_melhoresGenes])
geneNames=unlist(mget(geneFeatNames, illuminaHumanv4SYMBOL))
geneNames
geneFunction=unlist(mget(geneFeatNames, illuminaHumanv4GENENAME))
geneFunction

#######Controlo x Litio (baseline)#########
test2=rowttests(normalized_exp_filtrado2[normalized_exp_filtrado2$time=="baseline"],"agent")
test2
pvalues2= test2$p.value
pvalues2
pvalues2Signif= order(test2[test2$p.value<0.01, "p.value"])
pvalues2Signif
indices_melhoresGenes2=pvalues2Signif[1:25]
indices_melhoresGenes2
geneFeatNames2=featureNames(normalized_exp_filtrado2[indices_melhoresGenes2])
geneNames2=unlist(mget(geneFeatNames2, illuminaHumanv4SYMBOL))
geneNames2
geneFunction2=unlist(mget(geneFeatNames2, illuminaHumanv4GENENAME))
geneFunction2
>>>>>>> 0e15e90cf70421776e1fa0ec23f5e0eb4a6c4f07

#####Responde x Nao responde (baseline e litio)########
test3=rowttests(normalized_exp_filtrado2[normalized_exp_filtrado2$time=="baseline" & normalized_exp_filtrado2$agent=="lithium"],"other")
test3
pvalues3= test3$p.value
pvalues3
pvalues3Signif= order(test3[test3$p.value<0.01, "p.value"])
pvalues3Signif
indices_melhoresGenes3=pvalues3Signif[1:25]
indices_melhoresGenes3
geneFeatNames3=featureNames(normalized_exp_filtrado2[indices_melhoresGenes3])
geneNames3=unlist(mget(geneFeatNames3, illuminaHumanv4SYMBOL))
geneNames3
geneFunction3=unlist(mget(geneFeatNames3, illuminaHumanv4GENENAME))
geneFunction3

#####Responde x Nao responde (1month e litio)########
test4=rowttests(normalized_exp_filtrado2[normalized_exp_filtrado2$time=="1 month" & normalized_exp_filtrado2$agent=="lithium"],"other")
test4
pvalues4= test4$p.value
pvalues4
pvalues4Signif= order(test4[test4$p.value<0.01, "p.value"])
pvalues4Signif
indices_melhoresGenes4=pvalues4Signif[1:25]
indices_melhoresGenes4
geneFeatNames4=featureNames(normalized_exp_filtrado2[indices_melhoresGenes4])
geneNames4=unlist(mget(geneFeatNames4, illuminaHumanv4SYMBOL))
geneNames4
geneFunction4=unlist(mget(geneFeatNames4, illuminaHumanv4GENENAME))
geneFunction4

######Baseline x 1 month (Litio)##### 

test5=rowttests(normalized_exp_filtrado2[normalized_exp_filtrado2$agent=="lithium"],"time")
test5
pvalues5= test5$p.value
pvalues5
pvalues5Signif= order(test5[test5$p.value<0.05, "p.value"])
pvalues5Signif
indices_melhoresGenes5=pvalues5Signif[1:25]
indices_melhoresGenes5
geneFeatNames5=featureNames(normalized_exp_filtrado2[indices_melhoresGenes5])
geneNames5=unlist(mget(geneFeatNames5, illuminaHumanv4SYMBOL))
geneNames5
geneFunction5=unlist(mget(geneFeatNames5, illuminaHumanv4GENENAME))
geneFunction5


#####Baseline x 1 month (controlo)#####

test6=rowttests(normalized_exp_filtrado2[normalized_exp_filtrado2$agent=="control"],"time")
test6
pvalues6= test6$p.value
pvalues6
pvalues6Signif= order(test6[test6$p.value<0.05, "p.value"])
pvalues6Signif
indices_melhoresGenes6=pvalues6Signif[1:25]
indices_melhoresGenes6
geneFeatNames6=featureNames(normalized_exp_filtrado2[indices_melhoresGenes6])
geneNames6=unlist(mget(geneFeatNames6, illuminaHumanv4SYMBOL))
geneNames6
geneFunction6=unlist(mget(geneFeatNames6, illuminaHumanv4GENENAME))
geneFunction6

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
remove(gds2)
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
data$description[1:5]# [1] "Value for GSM1105438: OPT_1-baseline; src: Peripheral blood" 
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

#Histograma filtrado
sds__1=rowSds(exp_filtrado)
hist(sds__1,breaks=50,col='mistyrose')

####### Analise de expressao diferencial 
library('genefilter')
exp_filtrado2=data[sds>=2*m,]
test=rowttests(exp_filtrado2,"agent")
p_value=test$p.value#p_values do teste feito anteriomente
n_signifi=p_value[which(p_value<0.05)]#628 têm diferenças significativas

rank=order(test$p.value)#ordenar os p-values
ind_p25=rank[1:25]#?ndices dos 25 genes com menor p_value(maior evidencia de express?o diferencial)






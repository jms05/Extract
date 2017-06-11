########## Relatório ############ 

#******************carregar dados e visualizar os dados existentes**************************************
#install.packages("gplots")
#install.packages("dendextend")#para por por no heatmap nas linhas
#install.packages("dendextendRcpp")#para por por no heatmap nas linhas
#install.packages("d3heatmap")#para fazer o heatmap
source("http://bioconductor.org/biocLite.R")
biocLite()
#biocLite("GEOquery")
#biocLite("limma")
#biocLite("genefilter")

#biocLite("hgu133a.db")
library("hgu133a.db")
#biocLite("hgu95av2.db")

#biocLite("GOstats")
library(GOstats)
#biocLite("illuminaHumanv4.db")

library("illuminaHumanv4.db") #Get this library if you don't have

#biocLite("hgu133plus2.db")
#library("hgu133plus2.db")
#install.packages("stringr")

#install.packages("ggplot2")
#install.packages("ggdendro")

library(ggplot2)
library(ggdendro)

library(stringr)
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
#biocLite("affydata")
#biocLite("affy")
biocLite("affyPLM")
library(affydata)
library(affy)
library(affyPLM)
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

#####Baseline x 1 mês de tratamento (que respondem ao tratamento com lítio)########
#Na seguinte tabela encontram-se os genes com maior evidência estatística de expressão génica
#diferencial no que diz respeito ao atributo "time" nas condições de resposta e tratamento com lítio. 
#Este teste foi feito com o intuito de perceber, dentro dos indivíduos que foram tratados com lítio e responderam ao tratamento,
#quais os genes que após 1 mês do tratamento em questão (com lítio) apresentavam diferenças de expressão estatisticamente significativas
#relativamente à condição inicial (baseline). 
#De facto, constatou-se que respectivos p-values (p-values < 0.01) permitem rejeitar a hipótese nula, pelo
#que podemos assumir que os genes dos indivíduos que foram tratados com lítio e responderam ao tratamento são diferencialmente expressos
#no inicío do tratamento e após um mês do começo do mesmo.    

test3=rowttests(normalized_exp_filtrado2[normalized_exp_filtrado2$other=="responder" & normalized_exp_filtrado2$agent=="lithium"],"time")
test3
sign_pvalue3=test3[test3$p.value<0.01,]
sign_pvalue3
pvalues3Signif=sign_pvalue3[order(sign_pvalue3$p.value),]
nomes3=rownames(pvalues3Signif)
geneNames3=unlist(mget(nomes3, illuminaHumanv4SYMBOL))
geneNames3
geneFunction3=unlist(mget(nomes3, illuminaHumanv4GENENAME))
geneFunction3

#####Responde x Nao responde (1month e litio)########
#Na seguinte tabela encontram-se os 25 genes com maior evidência estatística de expressão génica
#diferencial no que diz respeito ao atributo "other" após 1 mês de tratamento com lítio. 
#Este teste foi feito com o intuito de perceber, dentro dos indivíduos que foram tratados com lítio durante 1 mês,
#quais os genes que apresentavam as maiores diferenças de expressão quando comparados os indivíduos que respondem e não respondem ao tratamento. 


test4=rowttests(normalized_exp_filtrado2[normalized_exp_filtrado2$time=="1 month" & normalized_exp_filtrado2$agent=="lithium"],"other")
test4
sign_pvalue4=test4[test4$p.value<0.01,]
sign_pvalue4
pvalues4Signif=sign_pvalue4[order(sign_pvalue4$p.value),]
nomes4=rownames(pvalues4Signif)[1:25]
geneNames4=unlist(mget(nomes4, illuminaHumanv4SYMBOL))
geneNames4
geneFunction4=unlist(mget(nomes4, illuminaHumanv4GENENAME))
geneFunction4

#####T-test para usar nos clustering
test5=rowttests(normalized_exp_filtrado2,"agent")
test5
sign_pvalue5=test5$p.value<0.05
sign_pvalue5

######### Analise de Enriquecimento 

#Funcao para ir buscar o ID do GO e o GO (gene ontology)
processaMarixRegEx <- function(data)
{
  matrix=as.matrix(data)
  names = c()
  ontology = c()
  index = 1
  regEX = "GO:[0-9]+"
  regexOnto = ".Ontology"
  nomes = row.names(matrix)
  for(i in 1:length(nomes)){
    val_regex = str_extract(nomes[i],regEX)
    if(!is.na(val_regex)){
      val_ontology = str_extract(nomes[i],regexOnto)
      if(!is.na(val_ontology)){
        names[index] = (val_regex)
        ontology[index] = (matrix[nomes[i],1])
        index = index+1
      }
    }
  }
  new.df= data.frame(GO = names,Ontology = ontology)
  return(new.df)
}

#Funcao para as contagens  de termos GO 
processTable <- function(tabela,data)
{
  chaves=as.character(data)
  names = c()
  contadores = c()
  for(i in 1:length(chaves)){
    contador = tabela[chaves[i]]
    names[i]=chaves[i]
    contadores[i] = contador
    
  }
  new.df= data.frame(
    GO = names,
    contadores = contadores)
  return(new.df)
}


#####Para o ttest1 - Controlo x Lítio (1 mês de tratamento) - "agent"

#2º titulo para o relatorio: 
#Analise de enriquecimento para os genes diferencialmente expressos 
#Após a realização do t-teste, procedemos à análise de enriquecimento a fim de obtermos os termos GO associados aos genes mais diferencialmente expressos.

filtrados1 = featureNames(normalized_exp_filtrado2) 
entrezID1= unlist(mget(filtrados1, illuminaHumanv4ENTREZID))  
pvalue1= test1$p.value<0.01
pvalFilt = normalized_exp_filtrado2[pvalue1,]
entrezIdSel1 = unlist(mget(featureNames(pvalFilt), illuminaHumanv4ENTREZID))
GO1 = new("GOHyperGParams", geneIds=entrezIdSel1, universeGeneIds=entrezID1, annotation="illuminaHumanv4.db", ontology="BP", pvalueCutoff=0.025, testDirection="over")
GO1_= hyperGTest(GO1)#Lista de grupos com menores p-values (ordem crescente)
summary(GO1_)

#Legenda: 
#Tabela com os 15 termos GO mais associados aos genes diferencialmente expressos no que diz respeito à distinção das amostras de controlo e lítio 
#ao fim de 1 mês de tratamento. 

#Na seguinte tabela, estão presentes os quinze termos GO com maior associação ao atributo "agent" (controlo/lítio), ao fim de 1 mês de tratamento. 
#O termo mais presente é "regulation of growth", o que nos indica que há genes, dentro dos mais diferencialmente expressos, envolvidos na regulação
#do crescimento celular. O lítio impede a formação de inositol, através da fosforilação de compostos, sendo responsável pelo aumento de Ca+ e Na+ no citoplasma neuronal, que por sua vez,
#nesta doença, aumenta o período refratário e a libertação de neurotransmissores, controlando assim, os sintomas. Assim, havendo um aumento destes iões 
#no citoplasma dos neurónios há, consequentemente, um aumento destes. Isto pode ser uma explicação para o facto de o termo GO mais presente, 
#dentro dos genes mais diferencialmente expressos, estar associado à regulação do crescimento celular, podendo isto ser um indício da influência 
#do lítio na expressão dos genes.  Livro Vera




#####Baseline x 1 mês de tratamento (que respondem ao tratamento com lítio)########
#2º titulo para o relatorio: 
#Analise de enriquecimento para os genes diferencialmente expressos
filtrados2 = featureNames(normalized_exp_filtrado2) 
entrezID2= unlist(mget(filtrados2, illuminaHumanv4ENTREZID))  
pvalue2= test3$p.value<0.01
pvalFilt = normalized_exp_filtrado2[pvalue2,]
entrezIdSel2 = unlist(mget(featureNames(pvalFilt), illuminaHumanv4ENTREZID))
GO2 = new("GOHyperGParams", geneIds=entrezIdSel2, universeGeneIds=entrezID2, annotation="illuminaHumanv4.db", ontology="BP", pvalueCutoff=0.025, testDirection="over")
GO2_= hyperGTest(GO2)#Lista de grupos com menores p-values (ordem crescente)
summary(GO2_)

#Legenda
#Tabela com os 15 termos GO mais associados aos genes diferencialmente expressos no que diz respeito à distinção das amostras baseline e primeiro mês
#de tratamento com lítio, em pacientes que responderam ao tratamento. 

#Na seguinte tabela, estão presentes os quinze termos GO com maior associação ao atributo "time" (baseline/1 mês de tratamento), em pacientes que responderam ao tratamento. 
#Sendo que nenhum dos termos se destacou no que diz respeito à contagem, não podemos tirar conclusões específicas sobre nenhum deles em particular. 
#No entanto, os termos resultantes estão relacionados com a fosforilação de compostos ("protein kinase A signaling") e associados a neurotransmissores
#("sequestering of neurotransmitter") e à regulação positiva da resposta humoral ("positive regulation of humoral immune response"). Estes resultados
#vão de encontro ao que é pretendido no tratamento da doença bipolar. Deste modo, podemos prever que estes genes, de pacientes que respondem ao tratamento com lítio,
#diferencialmente expressos nas condições de baseline e 1 mês de tratamento, estão relacionados com a resposta do organismo ao lítio.  Livro Vera

#####Responde x Nao responde (1month e litio)########
2º titulo para o relatorio: 
#Analise de enriquecimento para os genes diferencialmente expressos
filtrados3= featureNames(normalized_exp_filtrado2) 
entrezID3= unlist(mget(filtrados3, illuminaHumanv4ENTREZID))  
pvalue3= test4$p.value<0.01
pvalFilt = normalized_exp_filtrado2[pvalue3,]
entrezIdSel3 = unlist(mget(featureNames(pvalFilt), illuminaHumanv4ENTREZID))
GO3 = new("GOHyperGParams", geneIds=entrezIdSel3, universeGeneIds=entrezID3, annotation="illuminaHumanv4.db", ontology="BP", pvalueCutoff=0.025, testDirection="over")
GO3_= hyperGTest(GO3)#Lista de grupos com menores p-values (ordem crescente)
summary(GO3_)

#Legenda
#Tabela com os 15 termos GO mais associados aos genes diferencialmente expressos no que diz respeito à distinção das amostras dos pacientes que 
#respondem e não respondem ao tratamento no primeiro mês com lítio.

#Na seguinte tabela, estão presentes os quinze termos GO com maior associação ao atributo "other" (responde/não responde),
#em pacientes que foram tratados com lítio, ao fim de 1 mês.
#Os termos mais presentes, dentro dos genes mais diferencialmente expressos é "negative regulation of cell proliferation" e "negative regulation of cell migration".
#A descoberta da ocorrência de neurogénese (nascimento de novos neurónios) em cérebros adultos tem despoletado o interesse por estratégias de terapia 
#celular para o tratamento de doenças associadas ao sistema neuronal. A neurogénese, no cérebro adulto, ocorre em duas áreas: (a) Zona Subgranular 
#do hipocampo e (b) Zona Subventricular dos ventrículos laterais. 
#A partir desses locais, os novos neurónios migram em direção aos seus alvos finais em outras áreas cerebrais onde se diferenciam e 
#integram os circuitos locais. Diversos estudos têm mostrado que o tratamento farmacológico com lítio exerce efeitos sobre a neurogênese adulta, estimulando a sobrevivência e o 
#amadurecimento de novos neurónios gerados no cérebro de roedores adultos. 
#Assim, sendo que os genes com maior expressão diferencial, no que diz respeito a pacientes que responderam e não responderam ao tratamento com lítio,
#estão associados à proliferação e migração celular, podemos prever que o lítio seja capaz de despoletar estes mecanismos celulares e, deste modo, 
#controlar os sintomas da doença bipolar. 
#Artigo manhoso mas que da muito jeito 

#######CLUSTERING

#Funções utilizadas no Clustering 

#HeatMap
heatMapAnalysis <- function(data,colors)
{
  Rowv = data %>% dist %>% hclust %>% as.dendrogram %>%
    set("branches_k_color", k = 4) %>% set("branches_lwd", 4) %>% ladderize
  Colv = data %>% t %>% dist %>% hclust %>% as.dendrogram %>%
    set("branches_k_color", k = 2, value = colors) %>%
    set("branches_lwd", 4) %>% ladderize
  d3heatmap(data,  Colv = Colv, Rowv = Rowv)
}

#K-means 
k_means <- function(data_clustering,Ncenters,legends,lables=c("","")){
  #-> Kmeans  #Agent 
  dataFrame=exprs(data_clustering)
  km = kmeans(dataFrame, centers=Ncenters)  
  #gr?fico para explicar os resultados anteriores
  plot(dataFrame, col=km$cluster, pch=19, cex=2, xlab =lables[1], ylab=lables[2])
  points(km$centers, col=1:2, pch=3, cex=3, lwd=3)
  legend("topleft", legend = legends, cex = 1, pch = 19, col = 1:3)
  
}

#No que diz respeito à análise de dados, construiram-se dois dendrogramas, um para o 
#atributo "agent" e outro para o atributo "other", a fim de proceder ao estudo 
#dos clusters correspondentes. Para tal, usou-se apenas os vinte e cinco genes mais diferencialmente
#expressos do atributo dos respetivos metadados. A distância usada na obtenção
#dos clusters foi a euclideana. 

##Clustering associado ao agent 
test_agent=rowttests(normalized_exp_filtrado2,"agent")
sign_pvalue_agent=test_agent$p.value<0.01
data_clustering = normalized_exp_filtrado2[sign_pvalue_agent][1:25]
datacluster=exprs(data_clustering)
nomes_agent=rownames(normalized_exp_filtrado2[sign_pvalue_agent])[1:25]
clust_nomes_agent=unlist(mget(nomes_agent, illuminaHumanv4SYMBOL))

#Colocar os nomes dos genes no clustering 
rownames(datacluster) <- clust_nomes_agent

#Euclideana
dist_euc = dist(datacluster)
cl.hier.complete <- hclust(dist_euc)
plot(cl.hier.complete, xlab="Distância euclideana",ylab="Altura")

#Usando os vinte e cinco genes com a maior diferença de expressão, obtidos na análise
#de expressão diferencial deste atributo dos metadados, 
#é possível agrupar as amostras em dois grupos, tal como seria de prever, através da distância 
#entre os genes.

#Uma  representação  gráfica  dos  dados  relacionada  com  o  clustering  
#hierárquico  são os heatmaps.  
#Esta permite  representar  os  dados  como  uma  imagem  onde  cada  valor  da  
#matriz  de  dados  é  representado com uma célula cuja cor varia consoante o valor 
#respetivo, num gradiente  de cores que pode ser configurado. Os heatmaps tipicamente 
#incluem as árvores criadas  por  clustering  hierárquico  quer  ao  nível  das  
#linhas,  quer  ao  nível  das  colunas  de  dados. [Profs]

heatMapAnalysis(datacluster,c("orange", "blue")) #####Agent  

#A partir dos resultados do heatmap, foi possível verificar que os nossos genes
#se organizaram em dois clusters disitintos, tal como se observou no resultado do 
#clustering hierárquico. Assim, é possível deduzir que os genes relativos a cada 
#grupo, de controlo e de tratamento com lítio, apresentam, efetivamente, expressão 
#diferencial. 

k_means(data_clustering,2,c("Agent 1", "Agent 2")) ##Agent

#O  problema  de  clustering  k-means  constitui  uma  das  possíveis  formulações  
#do  clustering, onde o objetivo é o de minimizar a média do quadrado das distâncias 
#de cada ponto  ao  centro  do  cluster  a  que  pertence. 
#Assim, prevê-se que um dos grupos esteja mais intimamente ligado que outro. Isso pode também
#observar-se no clustering hierárquico, tendo em consideração as distâncias euclideanas.



##Clustering associado ao other
test_other=rowttests(normalized_exp_filtrado2,"other")
sign_pvalue_other=test_other$p.value<0.01
data_clustering2 = normalized_exp_filtrado2[sign_pvalue_other][1:25]
datacluster2=exprs(data_clustering2)
nomes_other=rownames(normalized_exp_filtrado2[sign_pvalue_other])[1:25]
clust_nomes_other=unlist(mget(nomes_other, illuminaHumanv4SYMBOL))

#Colocar os nomes dos genes no clustering 
rownames(datacluster2) <- clust_nomes_other

#Euclideana ##Melhor Resultado
dist_euc = dist(datacluster2)
cl.hier.complete <- hclust(dist_euc)
plot(cl.hier.complete, xlab="Distância euclideana",ylab="Altura")
#Usando os vinte e cinco genes com a maior diferença de expressão, obtidos na análise
#de expressão diferencial deste atributo dos metadados, 
#é possível agrupar as amostras em dois grupos, tal como seria de prever, através da distância 
#entre os genes.

heatMapAnalysis(datacluster2,c("orange", "blue")) ####Other 

#A partir dos resultados do heatmap, é possível verificar que a distribuição da 
#expressão dos genes de cada um dos grupos não é tão homogénea como se a que se observa no heatmap
#anterior. 
#Assim, isto pode indiciar que dentro dos genes dos indivíduos que respondem ao tratamento
#alguns deles podem expressar-se mais intensamente, participando, assim, mais de forma mais ativa
#na resposta ao tratamento. 

k_means(data_clustering2,2,c("Other 1", "Other 2")) ##Other 
#Tal como se verifica no heatmap anterior, a expressão génica dentro dos grupos é
#pouco homogénea o que sugere que nem todos os genes, que se expressam diferencialmente 
#no que diz respeito a condições de resposta ao tratamento, tenham a mesma
#participação na mesma, sendo previsível a existência de genes sub e sobrexpressos na resposta ao tratamento com lítio.




#######POR NO FINAL######
#Nos testes de hipóteses realizados concluímos que os genes com maiores diferenças estatísiticas significativas, obtidos para os t-test 
#realizado para os diferentes tipos de tratamento, não foram os mesmos para diferentes tempos experimentais. Isto sugere que, apesar 
#de não ser possível relacionar diretamente estas diferenças com o tratamento com lítio, diferentes tipos de tratamento podem estar relacionados
#com a variabilidade da expressão de certos genes. 
#Quanto às diferenças de expressão génica dentro de individuos que respondem ao tratamento com litio nas diferentes linhas temporais, verificamos
#que apenas 5 genes se expressavam de forma efetivamente diferente quando comparamos o tempo inicial do estudo com 1 mês de tratamento. 

#Relativamente ao teste estatístico que diz respeito à comparação de individuos que respondem e não respondem após 1 mês de tratamento 
#com lítio, os resultados insidem na possibilidade da resposta dos pacientes ao tratamento com lítio estar associada à regulação da expressão
#de alguns genes tornando-os, assim possíveis alvos para aumentar a eficácia do tratamento. 



#Instalação do GEOQuery
source("http://www.bioconductor.org/biocLite.R")
biocLite("GEOquery")

#Baixar o ficheiro GDS escolhido: GDS4766
library(Biobase)
library(GEOquery)
gds4766=getGEO('GDS4766', destdir=".")

#Baixar a anotação correspondente a este dataset:
biocLite("hgu133plus2.db")
library("hgu133plus2.db")


#PASSAR O OBJECTO GDS PARA UM OBJECTO EXPRESSION SET:
esetGDS4766=GDS2eSet(gds4766, do.log2=TRUE) #ao mesmo tempo que passa para um
            #objecto expression set, faz já um "passo" de pré-processamento:
            #transformação logarítmica
esetGDS4766
dim(esetGDS4766) #Dataset tem 54675 features (atributos) e 33 amostras
varMetadata(esetGDS4766) #Atributos do metadados das amostras:
                            #- sample (nome das amostras)
                            #- disease.state (se associado à gravidez ou não)
                            #- cell.type (epithelial ou stroma)
                            #- specimen (se normal ou tumor)
                            #- genotype (se receptor de estrogénio
                            #positivo ou negativo)
                            #- description (descrição de cada uma das amostras)


###PRÉ-PROCESSAMENTO###
#Não tem Valores omissos:
matExp=exprs(esetGDS4766)
sum(is.na(matExp))

#Flat Pattern Filters:
library(genefilter)
sds = rowSds(matExp)
m = median(sds)
hist(sds, breaks=50, col="mistyrose", main="Histograma dos desvios padrões",xlab="Desvio padrão",ylab="Frequência") 
abline(v=m, col="blue",lwd = 4,lty=3)
abline(v=m*2, col="red",lwd = 4,lty=3)
esetGDS4766.filtrado = esetGDS4766[sds >= 2*median(sds),]

#Normalização entre diferentes amostras
exp.filtrado=exprs(esetGDS4766.filtrado)
exprs(esetGDS4766.filtrado)=scale(exp.filtrado)





##########ANÁLISE DE EXPRESSÃO DIFERENCIAL##########
#identificar o conjunto de genes que têm níveis diferentes de expressão
#comparando duas condições experimentais
#Testes não-emparelhados: os dois grupos de amostras nunca se influenciam

#Análise entre os grupos associados à gravidez e não associados à gravidez
tt=rowttests(esetGDS4766.filtrado,"disease.state") 
rank.disease.state=order(tt$p.value)
p20.disease.state=rank.disease.state[1:20]
tt$p.value[p20.disease.state] #20 menores p.values
g.disease.state=featureNames(esetGDS4766.filtrado[p20.disease.state])
g.names.disease.state=unlist(mget(g.disease.state, hgu133plus2SYMBOL))

#Análise entre os grupos associados ao tipo de célula (estroma e epitelial)
tt2=rowttests(esetGDS4766.filtrado,"cell.type") 
rank.cell.type=order(tt2$p.value)
p20.cell.type=rank.cell.type[1:20]
tt2$p.value[p20.cell.type] #20 menores p.values
g.cell.type=featureNames(esetGDS4766.filtrado[p20.cell.type])
g.names.cell.type=unlist(mget(g.cell.type, hgu133plus2SYMBOL))

#Análise entre os grupos associados a specimen (tumor ou normal)
tt3=rowttests(esetGDS4766.filtrado,"specimen") 
rank.specimen=order(tt3$p.value)
p20.specimen=rank.specimen[1:20]
tt3$p.value[p20.specimen] #20 menores p.values
g.specimen=featureNames(esetGDS4766.filtrado[p20.specimen])
g.names.specimen=unlist(mget(g.specimen, hgu133plus2SYMBOL))

#Análise entre os grupos associados ao receptor de estrogénio (positivo ou negativo)
tt4=rowttests(esetGDS4766.filtrado,"genotype/variation") 
rank.genotype=order(tt4$p.value)
p20.genotype=rank.genotype[1:20]
tt4$p.value[p20.genotype] #20 menores p.values
g.genotype=featureNames(esetGDS4766.filtrado[p20.genotype])
g.names.genotype=unlist(mget(g.genotype, hgu133plus2SYMBOL))





########### Enrichment Analysis #############
library(GOstats)
filt = featureNames(esetGDS4766.filtrado) 
entrezUniverse= unlist(mget(filt, hgu133plus2ENTREZID)) 

#Análise entre os grupos associados à gravidez e não associados à gravidez
#Testes t para determinar genes diferencialmente expressos e seus IDs
pvalue = tt$p.value < 0.05 
pvalFilt = esetGDS4766.filtrado[pvalue,]
selectedEntrezIds = unlist(mget(featureNames(pvalFilt), hgu133plus2ENTREZID))
#Cria parâmetro e corre os testes estatísticos hipergeométricos:
#grupos de genes do Gene ontology, genes sobre expressos
params = new("GOHyperGParams", geneIds=selectedEntrezIds, universeGeneIds=entrezUniverse, annotation="hgu133plus2.db", ontology="MF", pvalueCutoff=0.025, testDirection="over")
hgOver = hyperGTest(params)#Lista de grupos com menores p-values (ordem crescente;
summary(hgOver)           #d? contagens de genes no grupo alvo e total de genes no grupo)

#An?lise entre os grupos associados ao tipo de c?lula (estroma e epitelial)
#Testes t para determinar genes diferencialmente expressos e seus IDs
pvalue2 = tt2$p.value < 0.05 
pvalFilt2 = esetGDS4766.filtrado[pvalue2,]
selectedEntrezIds2 = unlist(mget(featureNames(pvalFilt2), hgu133plus2ENTREZID))
#Cria parâmetro e corre os testes estatísticos hipergeométricos:
#grupos de genes do Gene ontology, genes sobre expressos
params2 = new("GOHyperGParams", geneIds=selectedEntrezIds2, universeGeneIds=entrezUniverse, annotation="hgu133plus2.db", ontology="MF", pvalueCutoff=0.025, testDirection="over")
hgOver2 = hyperGTest(params2)#Lista de grupos com menores p-values (ordem crescente;
summary(hgOver2)           #de contagens de genes no grupo alvo e total de genes no grupo)

#An?lise entre os grupos associados a specimen (tumor ou normal))
#Testes t para determinar genes diferencialmente expressos e seus IDs
pvalue3 = tt3$p.value < 0.05 
pvalFilt3 = esetGDS4766.filtrado[pvalue3,]
selectedEntrezIds3 = unlist(mget(featureNames(pvalFilt3), hgu133plus2ENTREZID))
#Cria parâmetro e corre os testes estatísticos hipergeométricos:
#grupos de genes do Gene ontology, genes sobre expressos
params3 = new("GOHyperGParams", geneIds=selectedEntrezIds3, universeGeneIds=entrezUniverse, annotation="hgu133plus2.db", ontology="MF", pvalueCutoff=0.025, testDirection="over")
hgOver3 = hyperGTest(params3)#Lista de grupos com menores p-values (ordem crescente;
summary(hgOver3)           #de contagens de genes no grupo alvo e total de genes no grupo)

#Análise entre os grupos associados ao receptor de estrog?nio (positivo ou negativo)
#Testes t para determinar genes diferencialmente expressos e seus IDs
pvalue4 = tt4$p.value < 0.05 
pvalFilt4 = esetGDS4766.filtrado[pvalue4,]
selectedEntrezIds4 = unlist(mget(featureNames(pvalFilt4), hgu133plus2ENTREZID))
#Cria parâmetro e corre os testes estatísticos hipergeométricos:
#grupos de genes do Gene ontology, genes sobre expressos
params4 = new("GOHyperGParams", geneIds=selectedEntrezIds4, universeGeneIds=entrezUniverse, annotation="hgu133plus2.db", ontology="MF", pvalueCutoff=0.025, testDirection="over")
hgOver4 = hyperGTest(params4)#Lista de grupos com menores p-values (ordem crescente;
summary(hgOver4)            #de contagens de genes no grupo alvo e total de genes no grupo)





#CLUSTERING

#Cluster com todas as amostras, para depois comparar se, de facto, escolher os top 20 genes
#ajuda a agrupar melhor as amostras.
euc.all=dist(t(exprs(esetGDS4766.filtrado)))
cl.hier.all=hclust(euc.all, method="average")

#Para os 20 genes mais "essenciais" para destinguir o grupo de gravidez/não gravidez
#Usando matriz euclediana e average
clust.p20.disease.state=esetGDS4766.filtrado[p20.disease.state]
euc.disease.state=dist(t(exprs(clust.p20.disease.state)))
cl.hier.disease.state=hclust(euc.disease.state, method="average")
#PLot com as cores segundo os dois grupos:
my.plot.hc(cl.hier.disease.state,lab.col=as.integer(esetGDS4766.filtrado$disease.state),lab=clust.p20.disease.state$sample,main="")
legend(21,5,legend=levels(esetGDS4766.filtrado$disease.state),text.col=c("black","red"),cex=0.6,bty="n",y.intersp=0.6,title="Legend:",title.adj=0,title.col="black")
#Correlação pearson e average
cor.disease.state=dist(1-cor(exprs(clust.p20.disease.state)))
cl.hier.cor.disease.state=hclust(cor.disease.state, method="average")
#PLot com as cores segundo os dois grupos:
my.plot.hc(cl.hier.cor.disease.state,lab.col=as.integer(esetGDS4766.filtrado$disease.state),lab=clust.p20.disease.state$sample)
legend(21,3.5,legend=levels(esetGDS4766.filtrado$disease.state),text.col=c("black","red"),cex=0.6,bty="n",y.intersp=0.6,title="Legend:",title.adj=0,title.col="black")
#Plot do cluster com todos os genes:
my.plot.hc(cl.hier.all,lab.col=as.integer(esetGDS4766.filtrado$disease.state),lab=esetGDS4766.filtrado$sample)
legend(21,90,legend=levels(esetGDS4766.filtrado$disease.state),text.col=c("black","red"),cex=0.6,bty="n",y.intersp=0.6,title="Legend:",title.adj=0,title.col="black")


#Para os 20 genes mais "essenciais" para destinguir o tipo de célula
#Usando matriz euclediana e average
clust.p20.cell.type=esetGDS4766.filtrado[p20.cell.type]
euc.cell.type=dist(t(exprs(clust.p20.cell.type)))
cl.hier.cell.type=hclust(euc.cell.type, method="average")
#PLot com as cores segundo os dois grupos:
my.plot.hc(cl.hier.cell.type,lab.col=as.integer(esetGDS4766.filtrado$cell.type),lab=clust.p20.cell.type$sample)
legend(25,5,legend=levels(esetGDS4766.filtrado$cell.type),text.col=c("black","red"),cex=0.6,bty="n",y.intersp=0.6,title="Legend:",title.adj=0,title.col="black")
#Correlação pearson e average
cor.cell.type=dist(1-cor(exprs(clust.p20.cell.type)))
cl.hier.cor.cell.type=hclust(cor.cell.type, method="average")
#PLot com as cores segundo os dois grupos:
my.plot.hc(cl.hier.cor.cell.type,lab.col=as.integer(esetGDS4766.filtrado$cell.type),lab=clust.p20.cell.type$sample)
legend(25,1.8,legend=levels(esetGDS4766.filtrado$cell.type),text.col=c("black","red"),cex=0.6,bty="n",y.intersp=0.6,title="Legend:",title.adj=0,title.col="black")
#Plot do cluster com todos os genes:
my.plot.hc(cl.hier.all,lab.col=as.integer(esetGDS4766.filtrado$cell.type),lab=esetGDS4766.filtrado$sample)
legend(25,90,legend=levels(esetGDS4766.filtrado$cell.type),text.col=c("black","red"),cex=0.6,bty="n",y.intersp=0.6,title="Legend:",title.adj=0,title.col="black")


#Para os 20 genes mais "essenciais" para destinguir o specimen
#Usando matriz euclediana e average
clust.p20.specimen=esetGDS4766.filtrado[p20.specimen]
euc.specimen=dist(t(exprs(clust.p20.specimen)))
cl.hier.specimen=hclust(euc.specimen, method="average")
#PLot com as cores segundo os dois grupos:
my.plot.hc(cl.hier.specimen,lab.col=as.integer(esetGDS4766.filtrado$specimen),lab=clust.p20.specimen$sample)
legend(25,7,legend=levels(esetGDS4766.filtrado$specimen),text.col=c("black","red"),cex=0.6,bty="n",y.intersp=0.6,title="Legend:",title.adj=0,title.col="black")
#Correlação pearson e average
cor.specimen=dist(1-cor(exprs(clust.p20.specimen)))
cl.hier.cor.specimen=hclust(cor.specimen, method="average")
#PLot com as cores segundo os dois grupos:
my.plot.hc(cl.hier.cor.specimen,lab.col=as.integer(esetGDS4766.filtrado$specimen),lab=clust.p20.specimen$sample)
legend(25,7,legend=levels(esetGDS4766.filtrado$specimen),text.col=c("black","red"),cex=0.6,bty="n",y.intersp=0.6,title="Legend:",title.adj=0,title.col="black")
#Plot do cluster com todos os genes:
my.plot.hc(cl.hier.all,lab.col=as.integer(esetGDS4766.filtrado$specimen),lab=esetGDS4766.filtrado$sample)
legend(25,90,legend=levels(esetGDS4766.filtrado$specimen),text.col=c("black","red"),cex=0.6,bty="n",y.intersp=0.6,title="Legend:",title.adj=0,title.col="black")


#Para os 20 genes mais "essenciais" para destinguir o receptor
#de estrogénio (positivo ou negativo)
#Usando matriz euclediana e average
clust.p20.genotype=esetGDS4766.filtrado[p20.genotype]
euc.genotype=dist(t(exprs(clust.p20.genotype)))
cl.hier.genotype=hclust(euc.genotype, method="average")
#PLot com as cores segundo os dois grupos:
my.plot.hc(cl.hier.genotype,lab.col=as.integer(esetGDS4766.filtrado$genotype),lab=clust.p20.genotype$sample)
legend(25,6,legend=levels(esetGDS4766.filtrado$genotype),text.col=c("black","red"),cex=0.6,bty="n",y.intersp=0.6,title="Legend:",title.adj=0,title.col="black")
#Correlação pearson e average
cor.genotype=dist(1-cor(exprs(clust.p20.genotype)))
cl.hier.cor.genotype=hclust(cor.genotype, method="average")
#PLot com as cores segundo os dois grupos:
my.plot.hc(cl.hier.cor.genotype,lab.col=as.integer(esetGDS4766.filtrado$genotype),lab=clust.p20.genotype$sample)
legend(25,4,legend=levels(esetGDS4766.filtrado$genotype),text.col=c("black","red"),cex=0.6,bty="n",y.intersp=0.6,title="Legend:",title.adj=0,title.col="black")
#Plot do cluster com todos os genes:
my.plot.hc(cl.hier.all,lab.col=as.integer(esetGDS4766.filtrado$genotype),lab=esetGDS4766.filtrado$sample)
legend(25,90,legend=levels(esetGDS4766.filtrado$genotype),text.col=c("black","red"),cex=0.6,bty="n",y.intersp=0.6,title="Legend:",title.adj=0,title.col="black")





########ANÁLISE PREDITIVA##########

#Vai-se fazer análise preditiva dos grupos disease.state e specimen

#Misturou-se a ordem das amostras, uma vez que em alguns atributos de factores
#(disease.state, por exemplo) têm as amostras por ordem ou seja, as amostras relativas
#a um dos níveis do factor aparecem primeiro e as amostras relativas a outro nível
#aparecem depois
set.seed(123)
aleatesetGDS4766.filtrado=esetGDS4766.filtrado[,sample(ncol(esetGDS4766.filtrado))]
#O número de exemplos de treino será 2/3 do número total de amostras
numExTreino=as.integer((2/3)*ncol(esetGDS4766.filtrado))


library(MLInterfaces)
##Modelos para prever a classe disease.state, se está associada à gravidez ou não
#MODELO K-Nearest Neighbours
knnResult.disease.state=MLearn(disease.state~., aleatesetGDS4766.filtrado, knnI(k=1), 1:numExTreino)
knn.confusionMat.disease.state=confuMat(knnResult.disease.state)
#A percentagem de exemplos correctamente classificados do modelo é 0.6363636:
accuracy.knnResult.disease.state=(knn.confusionMat.disease.state[1,1]+knn.confusionMat.disease.state[2,2])/sum(knn.confusionMat.disease.state)

#MODELO Naive-Bayes
nBayesResult.disease.state=MLearn(disease.state~., aleatesetGDS4766.filtrado, naiveBayesI, 1:numExTreino)
nBayes.confusionMat.disease.state=confuMat(nBayesResult.disease.state)
#A percentagem de exemplos correctamente classificados do modelo é 0.7272727:
accuracy.nBayesResult.disease.state=(nBayes.confusionMat.disease.state[1,1]+nBayes.confusionMat.disease.state[2,2])/sum(nBayes.confusionMat.disease.state)

#MODELO Árvore de decisão
arvoreResult.disease.state=MLearn(disease.state~., aleatesetGDS4766.filtrado, rpartI, 1:numExTreino)
arvore.confusionMat.disease.state=confuMat(arvoreResult.disease.state)
#A percentagem de exemplos correctamente classificados do modelo é 0.5454545:
accuracy.arvoreResult.disease.state=(arvore.confusionMat.disease.state[1,1]+arvore.confusionMat.disease.state[2,2])/sum(arvore.confusionMat.disease.state)
plot(RObject(arvoreResult.disease.state), uniform=T, branch=0.3,margin=0.1,compress=T)
text(RObject(arvoreResult.disease.state),use.n=T,cex=0.6)

#MODELO Random Forests (usa várias árvores, neste caso 600, de decisão para formar o modelo)
randomForResult.disease.state=MLearn(disease.state~., aleatesetGDS4766.filtrado, randomForestI, 1:numExTreino, ntree=600)
randomFor.confusionMat.disease.state=confuMat(randomForResult.disease.state)
#A percentagem de exemplos correctamente classificados do modelo é 0.7272727:
accuracy.randomForResult.disease.state=(randomFor.confusionMat.disease.state[1,1]+randomFor.confusionMat.disease.state[2,2])/sum(randomFor.confusionMat.disease.state)


###O PACKAGE MLInterfaces NÃO FAZ OPTIMIZAÇÃO DE PARÂMETROS, MAS O CARET FAZ:
library(caret)

data.mach.learn=data.frame(t(exprs(aleatesetGDS4766.filtrado))) #meter os dados num
                                          #data.frame com a matriz transpostas, para que
                                          #as amostras fiquem nas linhas e os atributos
                                          #nas colunas

#(os valores das accuracies podem mudar porque a forma como os grupos são escolhidos em
#crossed validation são aleatórios, daí fazer-se várias repetições (10), para que seja
#tida em conta a média das accuracies e não de apenas uma)

#MODELO knn, com optimização do parametro k (5 valores possíveis) e validação cruzada com
#fold 3. Todas as amostras foram, portanto, utilizadas para treinar e prever
knn.disease.state=train(data.mach.learn, aleatesetGDS4766.filtrado$disease.state, method="knn",tuneLength=5,trControl=trainControl("repeatedcv",number=3,repeats=10))
knn.disease.state
#O melhor k é 13, sendo este o valor final usado para o modelo e o accuracy médio 
#obtido é 0.597.
confusionMatrix(knn.disease.state)
#Variáveis mais importantes no modelo
varImp(knn.disease.state)
#nomes destes 20 genes mais importantes:
allMostImp.knn.disease.state=varImp(knn.disease.state)
top20Imp.knn.disease.state=na.exclude(data.frame(rownames(allMostImp.knn.disease.state$importance)[order(allMostImp.knn.disease.state$importance, decreasing=TRUE)]))[1:20,]
top20MostImp.knn.disease.state=mapply(substring,top20Imp.knn.disease.state,2, USE.NAMES = FALSE)
top20MostImp.knn.disease.state=unlist(mget(top20MostImp.knn.disease.state, hgu133plus2SYMBOL))


#MODELO naiveBayes, com validação cruzada com fold 3. Todas as amostras foram, portanto,
#utilizadas para treinar e prever
nBayes.disease.state=train(data.mach.learn, aleatesetGDS4766.filtrado$disease.state,method="nb",trControl=trainControl("repeatedcv",number=3,repeats=10))
nBayes.disease.state
#A accuracy média do modelo final é 0.6576, com o valor de usekernel falso
confusionMatrix(nBayes.disease.state)
#Variáveis mais importantes no modelo
allMostImp.nB=varImp(nBayes.disease.state)
allMostImp.nB
#Nomes dos genes das 20 variáveis mais importantes
top20Imp.nB.disease.state=na.exclude(data.frame(rownames(allMostImp.nB$importance)[order(allMostImp.nB$importance, decreasing=TRUE)]))[1:20,]
top20MostImp.nB.disease.state=mapply(substring,top20Imp.nB.disease.state,2, USE.NAMES = FALSE)
top20MostImp.nB.disease.state=unlist(mget(top20MostImp.nB.disease.state, hgu133plus2SYMBOL))

#MODELO Árvore de decisão, com optimização do parâmetro de complexidade (cp) (5 valores
#possíveis) e validação cruzada com fold 3. Todas as amostras foram, portanto, utilizadas
#para treinar e prever
arvore.disease.state=train(data.mach.learn, aleatesetGDS4766.filtrado$disease.state,method="rpart",tuneLength=5,trControl=trainControl("repeatedcv",number=3,repeats=10))
arvore.disease.state
#O melhor cp é 0.6923077, sendo este o valor final usado para o modelo e o accuracy
#médio é 0.6242.
#O modelo final obtido pode não ser um bom modelo, uma vez que a árvore de decisão apenas
#tem a raiz.
confusionMatrix(arvore.disease.state)
plot(arvore.disease.state) #A accuracy mantém-se com o aumento da complexidade
#Variáveis mais importantes no treino do modelo:
allMost.arv.disease.state=varImp(arvore.disease.state)

#MODELO random Forests, com optimização do parâmetro de complexidade (cp) (5 valores
#possíveis) e validação cruzada com fold 3. Todas as amostras foram, portanto, utilizadas
#para treinar e prever
randomForest.disease.state=train(data.mach.learn, aleatesetGDS4766.filtrado$disease.state,method="rf",trControl=trainControl("repeatedcv",number=3,repeats=10))
randomForest.disease.state
#O melhor modelo tem mtry=6378, sendo a accuracy médio do modelo final 0.633.
confusionMatrix(randomForest.disease.state)
plot(randomForest.disease.state) #A accuracy aumenta com o o número de atributos
                            #(seleccionados aleatoriamente)
#Variáveis mais importantes no treino do modelo:
allMostImp.forest.disease.state=varImp(randomForest.disease.state)
#Nomes dos genes das 20 variáveis mais importantes
top20Imp.forest.disease.state=na.exclude(data.frame(rownames(allMostImp.forest.disease.state$importance)[order(allMostImp.forest.disease.state$importance, decreasing=TRUE)]))[1:20,]
top20MostImp.forest.disease.state=mapply(substring,top20Imp.forest.disease.state,2, USE.NAMES = FALSE)
top20MostImp.forest.disease.state=unlist(mget(top20MostImp.forest.disease.state, hgu133plus2SYMBOL))

#MODELO pls, com optimização do parâmetro de número de componentes (5 valores possíveis
#e validação cruzada com fold 3. Todas as amostras foram, portanto, utilizadas para
#treinar e prever
pls.disease.state=train(data.mach.learn, aleatesetGDS4766.filtrado$disease.state,method="pls",tuneLength=5,trControl=trainControl("repeatedcv",number=3,repeats=10))
pls.disease.state
#O melhor parâmetro é 4, sendo este o valor final usado para o modelo e o accuracy médio
#é 0.7424.
confusionMatrix(pls.disease.state)
#20 Variáveis mais importantes no treino do modelo:
allMostImp.pls.disease.state=varImp(pls.disease.state)
#Nomes dos genes das 20 variáveis mais importantes
top20Imp.pls.disease.state=na.exclude(data.frame(rownames(allMostImp.pls.disease.state$importance)[order(allMostImp.pls.disease.state$importance, decreasing=TRUE)]))[1:20,]
top20MostImp.pls.disease.state=mapply(substring,top20Imp.pls.disease.state,2, USE.NAMES = FALSE)
top20MostImp.pls.disease.state=unlist(mget(top20MostImp.pls.disease.state, hgu133plus2SYMBOL))

#MODELO svm, com validação cruzada com fold 3. Todas as amostras foram portanto,
#utilizadas para treinar e prever
svm.disease.state=train(data.mach.learn, aleatesetGDS4766.filtrado$disease.state,method="svmLinear",trControl=trainControl("repeatedcv",number=3,repeats=10))
svm.disease.state
#A accuracy média do modelo final é  0.7394.
confusionMatrix(svm.disease.state)
#Variáveis mais importantes no treino do modelo:
allMostImp.svm.disease.state=varImp(svm.disease.state)
#Nomes dos genes das 20 variáveis mais importantes
top20Imp.svm.disease.state=na.exclude(data.frame(rownames(allMostImp.svm.disease.state$importance)[order(allMostImp.svm.disease.state$importance, decreasing=TRUE)]))[1:20,]
top20MostImp.svm.disease.state=mapply(substring,top20Imp.svm.disease.state,2, USE.NAMES = FALSE)
top20MostImp.svm.disease.state=unlist(mget(top20MostImp.svm.disease.state, hgu133plus2SYMBOL))

#Através das accuracies observadas, o modelo pls, com 4 componentes, parece ser o modelo
#que melhor se adequa a estimar se as amostras estão associadas à gravidez ou não.

####machine learning para tumor ou normal

#MODELO knn, com optimização do parametro k (5 valores possíveis) e validação cruzada com
#fold 3. Todas as amostras foram, portanto, utilizadas para treinar e prever
knn.specimen=train(data.mach.learn, aleatesetGDS4766.filtrado$specimen, method="knn",tuneLength=5,trControl=trainControl("repeatedcv",number=3,repeats=10))
knn.specimen
#O melhor k é 5, sendo este o valor final usado para o modelo e o accuracy médio 
#obtido é 0.9424.
confusionMatrix(knn.specimen)
#Variáveis mais importantes no modelo
allMostImp.knn.specimen=varImp(knn.specimen)
#Nomes dos genes das 20 variáveis mais importantes
top20Imp.knn.specimen=na.exclude(data.frame(rownames(allMostImp.knn.specimen$importance)[order(allMostImp.knn.specimen$importance, decreasing=TRUE)]))[1:20,]
top20MostImp.knn.specimen=mapply(substring,top20Imp.knn.specimen,2, USE.NAMES = FALSE)
top20MostImp.knn.specimen=unlist(mget(top20MostImp.knn.specimen, hgu133plus2SYMBOL))



#MODELO naiveBayes, com validação cruzada com fold 3. Todas as amostras foram, portanto,
#utilizadas para treinar e prever
nBayes.specimen=train(data.mach.learn, aleatesetGDS4766.filtrado$specimen,method="nb",trControl=trainControl("repeatedcv",number=3,repeats=10))
nBayes.specimen
#A accuracy média do modelo final é 0.9394, com o valor de usekernel falso
confusionMatrix(nBayes.specimen)
#Variáveis mais importantes no modelo
allMostImp.nB.specimen=varImp(nBayes.specimen)
allMostImp.nB.specimen
#Nomes dos genes das 20 variáveis mais importantes
top20Imp.nB.specimen=na.exclude(data.frame(rownames(allMostImp.nB.specimen$importance)[order(allMostImp.nB.specimen$importance, decreasing=TRUE)]))[1:20,]
top20MostImp.nB.specimen=mapply(substring,top20Imp.nB.specimen,2, USE.NAMES = FALSE)
top20MostImp.nB.specimen=unlist(mget(top20MostImp.nB.specimen, hgu133plus2SYMBOL))

#MODELO Árvore de decisão, com optimização do parâmetro de complexidade (cp) (5 valores
#possíveis) e validação cruzada com fold 3. Todas as amostras foram, portanto, utilizadas
#para treinar e prever
arvore.specimen=train(data.mach.learn, aleatesetGDS4766.filtrado$specimen,method="rpart",tuneLength=5,trControl=trainControl("repeatedcv",number=3,repeats=10))
arvore.specimen
#O melhor cp é 0.75, sendo este o valor final usado para o modelo e o accuracy
#médio é 0.9636.
#O modelo final obtido pode não ser um bom modelo, uma vez que a árvore de decisão apenas
#tem a raiz.
confusionMatrix(arvore.specimen)
plot(arvore.specimen) #A accuracy mantém-se até atingir uma complexidade perto de 0.8, a
            #partir da qual desce
#Plot da árvore:
plot(arvore.specimen$finalModel, uniform=T, branch=0.4,margin=0.1,compress=T)
text(arvore.specimen$finalModel, use.n=T,cex=0.9)
#Variáveis mais importantes no treino do modelo:
allMostImp.arv.specimen=varImp(arvore.specimen)
#Nomes dos genes das 20 variáveis mais importantes
top20Imp.arv.specimen=na.exclude(data.frame(rownames(allMostImp.arv.specimen$importance)[order(allMostImp.arv.specimen$importance, decreasing=TRUE)]))[1:20,]
top20MostImp.arv.specimen=mapply(substring,top20Imp.arv.specimen,2, USE.NAMES = FALSE)
top20MostImp.arv.specimen=unlist(mget(top20MostImp.arv.specimen, hgu133plus2SYMBOL))

#MODELO random Forests, com optimização do parâmetro de complexidade (cp) (5 valores
#possíveis) e validação cruzada com fold 3. Todas as amostras foram, portanto, utilizadas
#para treinar e prever
randomForest.specimen=train(data.mach.learn, aleatesetGDS4766.filtrado$specimen,method="rf",trControl=trainControl("repeatedcv",number=3,repeats=10))
randomForest.specimen
#O melhor modelo tem mtry=6378, sendo a accuracy médio do modelo final 0.9818.
confusionMatrix(randomForest.specimen)
plot(randomForest.specimen) #A accuracy aumenta com o o número de atributos
#(seleccionados aleatoriamente)
#Variáveis mais importantes no treino do modelo:
allMostImp.rf.specimen=varImp(randomForest.specimen)
#Nomes dos genes das 20 variáveis mais importantes
top20Imp.rf.specimen=na.exclude(data.frame(rownames(allMostImp.rf.specimen$importance)[order(allMostImp.rf.specimen$importance, decreasing=TRUE)]))[1:20,]
top20MostImp.rf.specimen=mapply(substring,top20Imp.rf.specimen,2, USE.NAMES = FALSE)
top20MostImp.rf.specimen=unlist(mget(top20MostImp.rf.specimen, hgu133plus2SYMBOL))

#MODELO pls, com optimização do parâmetro de número de componentes (5 valores possíveis
#e validação cruzada com fold 3. Todas as amostras foram, portanto, utilizadas para
#treinar e prever
pls.specimen=train(data.mach.learn, aleatesetGDS4766.filtrado$specimen,method="pls",tuneLength=10,trControl=trainControl("repeatedcv",number=3,repeats=10))
pls.specimen
#O melhor parâmetro é 5, sendo este o valor final usado para o modelo e o accuracy médio
#é 0.9879.
confusionMatrix(pls.specimen)
#20 Variáveis mais importantes no treino do modelo:
allMostImp.pls.specimen=varImp(pls.specimen)
#Nomes dos genes das 20 variáveis mais importantes
top20Imp.pls.specimen=na.exclude(data.frame(rownames(allMostImp.pls.specimen$importance)[order(allMostImp.pls.specimen$importance, decreasing=TRUE)]))[1:20,]
top20MostImp.pls.specimen=mapply(substring,top20Imp.pls.specimen,2, USE.NAMES = FALSE)
top20MostImp.pls.specimen=unlist(mget(top20MostImp.pls.specimen, hgu133plus2SYMBOL))

#MODELO svm, com validação cruzada com fold 3. Todas as amostras foram portanto,
#utilizadas para treinar e prever
svm.specimen=train(data.mach.learn, aleatesetGDS4766.filtrado$specimen,method="svmLinear",trControl=trainControl("repeatedcv",number=3,repeats=10))
svm.specimen
#A accuracy média do modelo final é  0.9667.
confusionMatrix(svm.specimen)
#Variáveis mais importantes no treino do modelo:
allMostImp.svm.specimen=varImp(svm.specimen)
#Nomes dos genes das 20 variáveis mais importantes
top20Imp.svm.specimen=na.exclude(data.frame(rownames(allMostImp.svm.specimen$importance)[order(allMostImp.svm.specimen$importance, decreasing=TRUE)]))[1:20,]
top20MostImp.svm.specimen=mapply(substring,top20Imp.svm.specimen,2, USE.NAMES = FALSE)
top20MostImp.svm.specimen=unlist(mget(top20MostImp.svm.specimen, hgu133plus2SYMBOL))

#Através das accuracies observadas, o modelo pls, com 5 componentes, parece ser, mais uma
#vez, o modelo que melhor se adequa a estimar se as amostras estão associadas a um tumor
#ou a células normais.





#####Feature Selection#####

control=rfeControl(functions=rfFuncs,method="repeatedcv",number=3, repeats=10)
#sizes é o tamanho de cada subconjunto a testar
x=floor(log2(6378))

#Selecção de atributos no caso de prever disease.state
selectedVariables=rfe(data.mach.learn,aleatesetGDS4766.filtrado$disease.state,sizes=2^(1:x),rfeControl=control)
selectedVariables #è um conjunto de 8 variáveis que é melhor capaz de prever
plot(selectedVariables,type=c("g","o")) #Plot que mostra a accuracy que é alcançada com
                #cada conjunto de variáveis, podendo-se observar que é com um conjunto
                #pequeno de variáveis que se conseguirá prever melhor
plot(selectedVariables,type=c("g","o"),xlim=range(0,10)) #Num plot com um limite no eixo
                #dos xx pode-se ver que, de facto, esse número de váriáveis é 8
variables.disease.state=predictors(selectedVariables) #dá o nome das 8 variáveis
#comparar com o top20 da expressão diferencial de disease.state??
#Repetir machine learning para este conjunto de variáveis??
#Saber o nome dos genes correspondentes às variáveis:
genes.fs1=mapply(substring,variables.disease.state,2, USE.NAMES = FALSE)
genes.fs2=unlist(mget(genes.fs, hgu133plus2SYMBOL))

#Selecção de atributos no caso de prever specimen
selectedVariables.specimen=rfe(data.mach.learn,aleatesetGDS4766.filtrado$specimen,sizes=2^(1:x),rfeControl=control)
selectedVariables.specimen #O conunto que melhor prevê é o com todas as variáveis
plot(selectedVariables.specimen,type=c("g","o"))



#####Respetição dos modelos realizados para disease.state, agora com o conjunto de 8
#####variáveis obtido

conjunto=data.mach.learn[,which(colnames(data.mach.learn)%in%variables.disease.state)]
conjuntoaleat=aleatesetGDS4766.filtrado[which(rownames(aleatesetGDS4766.filtrado)%in%genes.fs1),]

#MODELO knn, com optimização do parametro k (5 valores possíveis) e validação cruzada com
#fold 3. Todas as amostras foram, portanto, utilizadas para treinar e prever
knn.disease.state2=train(conjunto, conjuntoaleat$disease.state, method="knn",tuneLength=5,trControl=trainControl("repeatedcv",number=3,repeats=10))
knn.disease.state2
#O melhor k é 5, sendo este o valor final usado para o modelo e o accuracy médio 
#obtido é 0.9212.
confusionMatrix(knn.disease.state2)

#MODELO naiveBayes, com validação cruzada com fold 3. Todas as amostras foram, portanto,
#utilizadas para treinar e prever
nBayes.disease.state2=train(conjunto, conjuntoaleat$disease.state,method="nb",trControl=trainControl("repeatedcv",number=3,repeats=10))
nBayes.disease.state2
#A accuracy média do modelo final é 0.9424, com o valor de usekernel verdadeiro
confusionMatrix(nBayes.disease.state2)

#MODELO Árvore de decisão, com optimização do parâmetro de complexidade (cp) (5 valores
#possíveis) e validação cruzada com fold 3. Todas as amostras foram, portanto, utilizadas
#para treinar e prever
arvore.disease.state2=train(conjunto, conjuntoaleat$disease.state,method="rpart",tuneLength=5,trControl=trainControl("repeatedcv",number=3,repeats=10))
arvore.disease.state2
#O melhor cp é 0.6987879, sendo este o valor final usado para o modelo e o accuracy
#médio é 0.7.
#O modelo final obtido pode não ser um bom modelo, uma vez que a árvore de decisão apenas
#tem a raiz.
confusionMatrix(arvore.disease.state2)
plot(arvore.disease.state2) #A accuracy mantém-se com o aumento da complexidade

#MODELO random Forests, com optimização do parâmetro de complexidade (cp) (5 valores
#possíveis) e validação cruzada com fold 3. Todas as amostras foram, portanto, utilizadas
#para treinar e prever
randomForest.disease.state2=train(conjunto, conjuntoaleat$disease.state,method="rf",trControl=trainControl("repeatedcv",number=3,repeats=10))
randomForest.disease.state2
#O melhor modelo tem mtry=2, sendo a accuracy médio do modelo final 0.9545.
confusionMatrix(randomForest.disease.state2)
plot(randomForest.disease.state2) #A accuracy diminui com o o número de atributos
#(seleccionados aleatoriamente)

#MODELO pls, com optimização do parâmetro de número de componentes (5 valores possíveis
#e validação cruzada com fold 3. Todas as amostras foram, portanto, utilizadas para
#treinar e prever
pls.disease.state2=train(conjunto, conjuntoaleat$disease.state,method="pls",tuneLength=5,trControl=trainControl("repeatedcv",number=3,repeats=10))
pls.disease.state2
#O melhor parâmetro é 3, sendo este o valor final usado para o modelo e o accuracy médio
#é 0.9455.
confusionMatrix(pls.disease.state2)

#MODELO svm, com validação cruzada com fold 3. Todas as amostras foram portanto,
#utilizadas para treinar e prever
svm.disease.state2=train(conjunto, conjuntoaleat$disease.state,method="svmLinear",trControl=trainControl("repeatedcv",number=3,repeats=10))
svm.disease.state2
#A accuracy média do modelo final é  0.9152.
confusionMatrix(svm.disease.state2)

#Os valores das accuracies para todos os modelos de disease.state estão, de facto, maiores
#e o tempo que decorre para treinar um modelo é bem menor.

#Através das accuracies observadas, o modelo random forests, com mty 2, parece ser
#o modelo que melhor se adequa a estimar se as amostras estão associadas à gravidez ou
#não, com os novos atributos seleccionados. É preciso ter em atenção que a função
#escolhida para fazer a selecção de atributos foi rfFuncs (random forests??).





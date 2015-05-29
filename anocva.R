library(cluster)

anocva<-function(A,Aj,l){ #funcao recebe as dissimilaridades (de cada populacao e total) e a clusterizacao e retorno as estatisticas do teste
  k = length(Aj) #numero de populacoes
  N = ncol(A) #N
  
  #Abaixo o vetor S das silhuetas.
  S =  silhouette(l,dmatrix = A)[,3]
  #Abaixo uma lista com os vetor Sj de cada silhueta das populacoes.
  Sj = lapply(Aj,fun<-function(x){silhouette(l,x)[,3]})
  
  #Abaixo o DeltaS - estatistica do teste que testa se ha ao menos uma diferenca
  DeltaS = 0
  for (j in 1:k){
    DeltaS = DeltaS + sum((S-Sj[[j]])*(S-Sj[[j]]))
  }  
  
  #Abaixo o deltaS - estatistica do teste que testa em qual regiao ha diferenca
  deltaS = c()
  for (q in 1:N){
    temp_deltaS = 0
    for (j in 1:k){
      temp_deltaS = temp_deltaS + Sj[[j]][q]
    }
    temp_deltaS = temp_deltaS/k
    deltaS = c(deltaS,(S[q] - temp_deltaS)^2)
  }
  return(list(DeltaS,deltaS))
}

boot_strap <- function(dados,straps,fun_cluster,fun_ncluster,fun_dist = dist){#funcao que replica a funcao anova acima para um conjunto de dados nos straps
  k = length(dados) #numero de populacoes
  n = c();for(j in 1:k){n = c(n,length(dados[[j]]))} #tamanho de cada populacao. vetor de tamanho k
  lista = list();for (i in 1:length(n)){for (j in 1:n[i]){lista = c(lista,list(c(i,j)))}} #crio uma lista com pares individuos/populacao para ajudar no boot strap
  
  Aij = lapply(dados,fun<-function(x){lapply(x,fun2<-function(x){as.matrix(fun_dist(t(x)))})}) #dissimilaridades de cada individuo. usa qualquer fun_dist para calucular a distancia
  Aj = lapply(Aij,fun<-function(x){X=0;for(i in x){X = X + i};return(X/length(x))}) #dissimilaridades de cada populacao
  A = 0;for(j in 1:k){A = A + Aj[[j]] * n[j] / sum(n)} #dissimilaridades de todo o conjunto
  c = fun_ncluster(A,fun_cluster)
  l = fun_cluster(A,c) #funcao que difene os clusters
  S = anocva(A,Aj,l) #calculo as estatisticas de teste para o conjunto de dados
  S_temp = list() #lista que vai carregar as estatisticas de teste para o boot straps
  for (strap in 1:straps){
    fun<-function(x){Aij[[lista[[x]][1]]][[lista[[x]][2]]]} #funcao que retorna as dissimilaridades dos straps
    fun2<-function(x){lapply(sample(1:sum(n),n[x],replace = TRUE),fun)} #funcao que seleciona aleatoriamente os straps
    Aij_boot = lapply((1:k),fun2) #Aij_boot carrega as dissimilaridades desse strap. Note que nao sorteio os dados, mas as dissimilaridades. Dessa maneira nao preciso calcular as dissimilaridades varias vezes, que eh a parte critica da funcao.
    Aj = lapply(Aij_boot,fun<-function(x){X=0;for(i in x){X = X + i};return(X/length(x))}) #mesmo que pra populacao
    A = 0;for(j in 1:k){A = A + Aj[[j]] * n[j] / sum(n)} #mesmo que pra populacao
    c = fun_ncluster(A,fun_cluster)
    l = fun_cluster(A,c) #mesmo que pra populacao
    S_temp = c(S_temp,list(anocva(A,Aj,l))) #mesmo que pra populacao
  }
  #O que a funcao retorna: uma lista com dois elementos. O primeiro eh o valor-p do teste que compara se ao menos uma populacao apresenta diferenca
  #o segundo valor da lista eh um vetor com N p-valores, um para cada regiao.
  return(list(mean(sapply(X = S_temp,FUN<-function(X) S[[1]]<X[[1]])),apply(sapply(X = S_temp,FUN<-function(X) S[[2]]<X[[2]]),MARGIN = 1,mean))) #comparo os straps com os dados originais.
}
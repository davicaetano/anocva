library(MASS)
amostra <- function(n1,delta,N1=10,N2=1,sigma1 = 1,sigma2=1){ #funcao que gera os dados para a simulacao descrita em 2.4
  Sigma1 = rbind(c(sigma1,0),c(0,sigma1))
  Sigma2 = rbind(c(sigma2,0),c(0,sigma2))
  Sigma3 = rbind(c(sigma1,0),c(0,sigma1))
  N3 = N1
  normal = list()
  adhd = list()
  for (i in 1:n1){
    p1 = t(rbind(mvrnorm(N1,c(0,0),Sigma1),mvrnorm(N2,c(0+delta*1,0),Sigma2),mvrnorm(N3,c(2,0),Sigma3)))
    normal = c(normal,list(p1))
    p2 = t(rbind(mvrnorm(N1,c(0,0),Sigma1),mvrnorm(N2,c(0+delta*0,0),Sigma2),mvrnorm(N3,c(2,0),Sigma3)))
    adhd = c(adhd,list(p2))
  } 
  return(list(adhd,normal))
}

#gera num_am monte carlo simulacoes para cada combinacao deltas x n1s
#Eh necessario passar como parametros as funcoes que clusteriza e a funcao que determina o numero de clusters
#A funcao que clusteriza deve retornar um vetor de tamanho N com os numeros dos clusters de cada regiao q
#A funcao que define o numero de clusters pode ser uma constante, mas nao eh aconselhavel que seja
simula<-function(deltas,n1s,straps,num_am,fun_cluster,fun_ncluster){
  p1 = list()
  straps = 300
  delta = deltas[1]
  for (delta in deltas){
    p2 = list()
    n1 = n1s[1]
    for (n1 in n1s){
      p3 = list()
      for(am in 1:num_am){
        dados = amostra(n1,delta)
        p3 = c(p3,list(boot_strap(dados,straps,fun_cluster,fun_ncluster)))
      }
      p2 = c(p2,list(p3))
    }
    p1 = c(p1,list(p2))
  }
  return(p1)
}
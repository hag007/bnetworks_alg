library("COSINE")


diff_gen_PPI <-
  function(data1,data2,PPI,pvals=NULL){
    num_sample_1 <- dim(data1)[1]
    num_sample_2 <- dim(data2)[1]
    num_gene <- dim(data1)[2]
    PPI <-PPI[PPI[,1]%in% colnames(data1) & PPI[,2]%in% colnames(data1),]
    num_edge <- dim(PPI)[1]
    node_score <- rep(0,num_gene)
    edge_score <- rep(0,num_edge)
    type <- c(rep(0,num_sample_1),rep(1,num_sample_2))
    
    # calculate the statistics measuing the differential expression of each gene between the 2 groups
    
    for(i in 1:num_gene){
      if (is.null(pvals)){
        data <- c(data1[,i],data2[,i])
        node_score[i] <- f.test(data,type)
      }
      else{
        node_score[i] <- -log10(pvals[i])
      }
      print(i)
    }
    node_score[is.na(node_score)] <- mean(node_score[!is.na(node_score)])

    
    # calculate the statistics measuing the differential co-expression of each gene-pair between the 2 groups
    
    for(i in 1:num_edge){
      gene1 <- as.character(PPI[i,1])
      gene2 <- as.character(PPI[i,2])
      data.x <- c(data1[,gene1],data2[,gene1])
      data.y <- c(data1[,gene2],data2[,gene2])
      edge_score[i] <- cov(data.x,data.y)
      print(i)
    }
    
    edge_score[(edge_score==0)] <- mean(edge_score[!(edge_score==0)])
    
    scaled_node_score <- (node_score - mean(node_score))/sd(node_score)
    scaled_edge_score <- (edge_score - mean(edge_score))/sd(edge_score)
    names(scaled_node_score) <- colnames(data1)
    return(list(scaled_node_score, scaled_edge_score, PPI))
  }



get_quantiles_PPI <-
  function(scaled_node_score,scaled_edge_score,PPI,klist,pop_size){
    
    
    # The function to get the "node_score_term" and "edge_score_term" of a sub-network denoted by "vector"
    # "vector" is a binary vector with length equal to the size of the whole network. 
    # An element of value "1" indicates the inclusion of that gene in the selected sub-network.
    
    
    node_edge<-function(sub){    
      n<-length(sub)    
      node_score<-sum(scaled_node_score[sub])/sqrt(n)
      edges<- PPI[,1] %in% sub & PPI[,2] %in% sub
      m<-sum(edges)
      edge_score<-sum(scaled_edge_score[edges])/sqrt(m)
      return(c(node_score,edge_score))
    }
    
    n<-length(scaled_node_score)
    all_genes<-names(scaled_node_score)
    edge_score_term<-vector(length=length(klist),mode="list")
    node_score_term<-vector(length=length(klist),mode="list")
    save_sub<-NULL
    for(i in 1:length(klist)){
      k<-klist[i]
      for(j in 1:pop_size){         
        sub<-random_network_sampling_PPI(k,PPI,all_genes)
        node_edge_score<-node_edge(sub)
        node_score_term[[i]][j]<-node_edge_score[1]
        edge_score_term[[i]][j]<-node_edge_score[2]
        if(is.na(edge_score_term[[i]][j])) save_sub<-sub
        print(c(i,j))
      }   
      
    }
    
    
    log_abs_edge_node_ratio <- vector(length=length(klist),mode="list")
    for(i in 1:length(klist)){
      log_abs_edge_node_ratio[[i]] <- log10(abs(edge_score_term[[i]]/node_score_term[[i]]))
    }
    b <- NULL
    for(i in 1:length(klist)){
      a <- summary(log_abs_edge_node_ratio[[i]])
      b <- rbind(b,a)
    }
    
    ratio <- apply(b,2,mean)[-4]
    lambda <- sort(1/(1+10^ratio))
    names(lambda)<-names(ratio)
    return(list(ratio,lambda))
  }





# data("scaled_node_score")
# data("scaled_edge_score")
# data(PPI)

# data(simulated_data)
# data1 <- simulated_data[[1]]
# data2 <- simulated_data[[7]]
# colnames(data1)<-colnames(data2)<-as.character(1:500)
# test <- diff_gen_PPI(data1, data2, PPI)
# 
# scaled_node_score <- test[[1]]
# scaled_edge_score <- test[[2]]




PPI <- as.matrix(read.table("/home/hag007/bnet/output/ppi_i.txt", stringsAsFactors = F, header=T))
ge <- as.matrix(read.table("/home/hag007/bnet/datasets/TNFa_2/output/cosine_ge.tsv", stringsAsFactors = F, row.names=1 ,header=T))
ge <-t(ge)
ge.control <- ge[c(1,2,3),]
ge.case <- ge[c(4,5,6),]

colnames(ge.case)<-colnames(ge.control)<-as.character(1:dim(ge.case)[2])
test <- diff_gen_PPI(ge.case,ge.control, PPI,pvals=ge[7,] ) #

scaled_node_score <- test[[1]]
scaled_edge_score <- test[[2]]



quantiles<-get_quantiles_PPI(scaled_node_score,scaled_edge_score, PPI, klist=seq(5,15,by=5),pop_size=10000)

q.ratio <- quantiles[[1]]
q.lambda <- quantiles[[2]]

bs = c()
ss = c()
adjs = c()
scores = c()
subnets = list(rep(NULL,5))

for (l in 1:5){

  GA_result<-GA_search_PPI(lambda=q.lambda[[l]],scaled_node_score,scaled_edge_score, PPI, 
                           num_iter=300, muCh=0.005, zToR=7, minsize=1) #
  bs <- c(bs, GA_result$Best_Scores)
  ss <- c(ss, GA_result$Subnet_size)
  subnets[[l]] <- list(0)
  subnets[[l]] <- GA_result$Subnet
  as = Score_adjust_PPI(scaled_node_score, scaled_edge_score,
                        PPI, q.lambda[[l]], GA_result$Subnet, 1000, GA_result$Best_Scores)
  scores <- c(scores, GA_result$Best_Scores)
  adjs <- c(adjs, as)
}

max.subnet = subnets[[match(max(adjs),adjs)]]

# chosen.lambda <- choose_lambda(diff_expr=scaled_node_score1,
#                                     diff_coex=scaled_edge_score1,q.lambda,subnet_size=ss,
#                                     num_random_sampling=10000,best_score=bs)





cmini<-function(data,lamda,ord=NULL){
  print(lamda)
  n_gene<-length(data[1,])
  G<-matrix(1,n_gene,n_gene)
  diag(G)<-0
  Gval<-G
  order<--1
  t<-0
  while(t==0){
    order <- order + 1
    if(!is.null(ord)){
      if(order>ord){
        order<-order - 1
        stop("Order Crossed!")
      }
    }
    a<-pathConsistency(G,Gval,order,data,t,lamda)
    G<-a$G
    Gval<-a$Gval
    t<-a$t
    #print(Gval)
    if(t==0)
    {
      break
      #stop("No edge is reduce")
    }else{t<-0}
  }
  order<-order - 1
  #print(Gval)
  return(a)
}

pathConsistency<-function(G,Gval,order,data,t,lamda){

  if(order==0){
    for(i in 1:length(G[1,])){
      for(j in 1:length(G[1,])){
        if(G[i,j]!=0){
          cmiv<-cmi(data[i,],data[j,]);
          Gval[i,j]<-cmiv
          Gval[j,i]<-cmiv
          if(cmiv<lamda){
            G[i,j]<-0
            G[j,i]<-0
          }
        }
      }
    }
    #print(G)
    t<-t+1
  }else{
    for(i in 1:length(G[1,])){
      for(j in 1:length(G[1,])){
        if(G[i,j]!=0){
          adj<-c()
          for(k in 1:length(G[1,])){
            if(G[i,k]!=0 && G[j,k]!=0){
              adj<-c(adj,k)
            }
          }

          if(length(adj)>=order){
            if(length(adj)==1){
              combntnslist<-combn(adj,order)
              combntnslist<-matrix(combntnslist[,length(combntnslist)],1,1)
            }else{
              combntnslist<-t(combn(adj,order))
            }
            combntnsrow<-length(combntnslist)/order
            cmiv<-0
            v1<-data[i,]
            v2<-data[j,]
            #print(combntnslist)
            for(k in 1:combntnsrow){
              vcs<-data[combntnslist[k,],]
              a<-MI2(v1,v2,vcs)
              cmiv<-max(cmiv,a)
            }
            #print(cmiv)
            Gval[i,j]<-cmiv
            Gval[j,i]<-cmiv
            if(cmiv<lamda){
              G[i,j]<-0
              G[j,i]<-0
            }
            t<-t + 1
          }
        }
      }
    }
  }
  retlst<-list("G"=G,"Gval"= Gval, "t"=t)
  return (retlst)
}
cmi<-function(v1,v2,vcs=NULL){
  if(is.null(vcs)){
    c1<-det(cov(matrix(v1)))
    c2<-det(cov(matrix(v2)))
    ct<-c(v1,v2)
    ct<-matrix(ct,length(v1),2)
    c3<-det(cov(ct))
    cmiv<-0.5*log(c1*c2/c3)
  }else if(!is.null(vcs)){
    #need to implement
  }
  if(is.nan(cmiv)){
    cmiv<-1.0e+10
  }
  return(cmiv)
}

MI2<-function(x,y,z){
  r_dmi<-(cas(x,y,z)+cas(y,x,z))/2;
  return(r_dmi)
}
cas<-function(x,y,z){
  #print(x)
  #print(y)
  #print(z)
  z<-matrix(z,1,length(z))
  n1<-length(z[,1])
  n<-n1+2
  Cov<-cov(matrix(x,byrow = TRUE))
  #print(Cov)
  xyz<-c(x,y,z)
  Covm<-cov(t(matrix(xyz,3,length(x),byrow = TRUE)))
  #print(Covm)
  xz<-c(x,z)
  Covm1<-cov(t(matrix(xz,2,length(x),byrow = TRUE)))
  #print(Covm1)

  InvCov<-solve(Cov)
  InvCovm<-solve(Covm)
  InvCovm1<-solve(Covm1)
  #print("Here at InCovm1")
  C11 <- InvCovm1[1,1]
  C12 <- 0
  #print(InvCovm1)
  C13 <- InvCovm1[1,2:(n1+1)]
  #print("Here at InCovm1")
  #print(InvCovm)
  #print(InvCov)
  #print(InvCov[1,1])
  C23 <- InvCovm[2,3:(2+n1)]-InvCovm[1,2] * (1/(InvCovm[1,1]-InvCovm1[1,1]+InvCov[1,1])) * (InvCovm[1,3:(2+n1)] - InvCovm1[1,2:(1+n1)])
  C22 <- InvCovm[2,2]- (InvCovm[1,2] %*% InvCovm[1,2]) %*% (1/(InvCovm[1,1]-InvCovm1[1,1]+InvCov[1,1]))
  C33 <- InvCovm[3:(2+n1),3:(2+n1)]- (1/(InvCovm[1,1]-InvCovm1[1,1]+InvCov[1,1])) %*% (t(InvCovm[1,3:(2+n1)]-InvCovm1[1,2:(1+n1)])%*%(InvCovm[1,3:(2+n1)]-InvCovm1[1,2:(1+n1)]))
  #print(C23)
  #print(C22)
  #print(C33)

  c1123 <- c(C11,C12,C13)
  c2123 <- c(C12,C22,C23)
  c3123 <- c(t(C13),t(C23),C33)
  ctmp <- c(c1123,c2123,c3123)
  InvC <- matrix(ctmp,3,length(c1123),byrow = TRUE)

  C0 <- Cov[1,1] * (InvCovm[1,1] - InvCovm1[1,1] + InvCov[1,1])
  CS <-0.5*(sum(diag((InvC %*% Covm)))+log(C0)-n)
  return(CS)
}

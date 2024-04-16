## Create the matrix A that converts between observations and DID estimators
### Several helper functions are used to create the components of A

gen_Aj <- function(j) {
  return(cbind(matrix(data=rep(-1, j-1), nrow=j-1, ncol=1),
               diag(nrow=j-1)))
}

gen_Adot <- function(J) {
  do.call(rbind,
          lapply(X=J:2,
                 FUN=function(x) cbind(matrix(data=0, nrow=x-1, ncol=J-x),
                                       gen_Aj(x))))
}

gen_Arow <- function(N,J,n1,n2,J2,Adot) {
  cbind(matrix(data=0, nrow=J2, ncol=(n1-1)*J),
        Adot,
        matrix(data=0, nrow=J2, ncol=(n2-n1-1)*J),
        Adot,
        matrix(data=0, nrow=J2, ncol=(N-n2)*J))
}

gen_A <- function(N,J) {
  J2 <- J*(J-1)/2
  Adot <- gen_Adot(J)
  
  do.call(rbind,
          lapply(X=1:(N-1),
                 FUN=function(n1) do.call(rbind,
                                          lapply(X=(n1+1):N,
                                                 FUN=function(n2) gen_Arow(N,J,n1,n2,J2,Adot)))))
}

## Create a data frame where each row corresponds to an entry in vector D
### The columns give i, i', j, and j' corresponding to that D
gen_D <- function(N,J) {
  Cl <- data.frame(i=unlist(sapply(1:(N-1), FUN=function(x) rep(x,N-x))),
                   i.prime=unlist(sapply(1:(N-1), FUN=function(x) (x+1):N)))
  Pd <- data.frame(j=unlist(sapply(1:(J-1), FUN=function(x) rep(x,J-x))),
                   j.prime=unlist(sapply(1:(J-1), FUN=function(x) (x+1):J)))
  return(cross_join(Cl,Pd))
}

## Old versions of the gen_A function:
# gen_A <- function(N,J) {
#   J2 <- J*(J-1)/2
#   Adot <- gen_Adot(J)
#   
#   A <- NULL
#   for (n1 in 1:(N-1)) {
#     for (n2 in (n1+1):N) {
#       A <- rbind(A, gen_Arow(N,J,n1,n2,J2,Adot))
#     }
#   }
#   return(A)
# }
# 
# gen_A2 <- function(N,J) {
#   J2 <- J*(J-1)/2
#   Adot <- gen_Adot(J)
#   
#   A <- NULL
#   for (n1 in 1:(N-1)) {
#     A <- rbind(A, do.call(rbind,
#                           lapply(X=(n1+1):N,
#                                  FUN=function(n2) gen_Arow(N,J,n1,n2,J2,Adot))))
#   }
#   return(A)
# }


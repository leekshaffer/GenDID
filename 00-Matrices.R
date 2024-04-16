### Functions to create A matrix ###
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

gen_Arow <- function(J,N,n1,n2,J2,Adot) {
  cbind(matrix(data=0, nrow=J2, ncol=(n1-1)*J),
        Adot,
        matrix(data=0, nrow=J2, ncol=(n2-n1-1)*J),
        Adot,
        matrix(data=0, nrow=J2, ncol=(N-n2)*J))
}

gen_A <- function(J,N) {
  J2 <- J*(J-1)/2
  Adot <- gen_Adot(J)
  
  do.call(rbind,
          lapply(X=1:(N-1),
                 FUN=function(n1) do.call(rbind,
                                          lapply(X=(n1+1):N,
                                                 FUN=function(n2) gen_Arow(J,N,n1,n2,J2,Adot)))))
}

## Old versions of the gen_A function:
# gen_A <- function(J,N) {
#   J2 <- J*(J-1)/2
#   Adot <- gen_Adot(J)
#   
#   A <- NULL
#   for (n1 in 1:(N-1)) {
#     for (n2 in (n1+1):N) {
#       A <- rbind(A, gen_Arow(J,N,n1,n2,J2,Adot))
#     }
#   }
#   return(A)
# }
# 
# gen_A2 <- function(J,N) {
#   J2 <- J*(J-1)/2
#   Adot <- gen_Adot(J)
#   
#   A <- NULL
#   for (n1 in 1:(N-1)) {
#     A <- rbind(A, do.call(rbind,
#                           lapply(X=(n1+1):N,
#                                  FUN=function(n2) gen_Arow(J,N,n1,n2,J2,Adot))))
#   }
#   return(A)
# }


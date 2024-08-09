#######################################
###### File: A_Const.R ################
###### Lee Kennedy-Shaffer ############
###### Created 2024/04/16 #############
###### Updated 2024/08/08 #############
#######################################

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
        -1*Adot,
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
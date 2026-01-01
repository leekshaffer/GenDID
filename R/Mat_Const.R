#######################################
###### File: Mat_Const.R ##############
###### Lee Kennedy-Shaffer ############
#######################################

# Load packages

require(dplyr)
require(tibble)
require(tidyr)

# gen_js function

## Main gen_js function

### Inputs:
#### Clusters: a vector of the unique unit labels
#### StartPeriods: a vector of the starting periods (period labels) corresponding
####  to each unit in Clusters, in the same order
#### OrderedPds: a vector of the temporal order of the period labels that will
####  have outcomes (may have more periods than StartPeriods)
### Output: A list of the following:
#### Start_js: DF of units ordered by their start period,
####  with numeric labels corresponding to their order
#### PdLabels: DF of the period labels and their numbered order
#### N: number of periods
#### J: number of units

gen_js <- function(Clusters,StartPeriods,OrderedPds) {
  if (length(unique(OrderedPds)) != length(OrderedPds)) {
    warning(simpleWarning("There are repeated names in OrderedPds. The first occurrence only will be used."))
    OrderedPds <- unique(OrderedPds)
  }
  J <- length(OrderedPds)

  if (length(Clusters) != length(unique(Clusters))) {
    stop(simpleError("Each cluster must appear only once."))
  }
  N <- length(unique(Clusters))
  if (N < 2) {
    stop(simpleError("There must be at least 2 clusters."))
  } else if (N != length(StartPeriods)) {
    stop(simpleError("The length of Clusters and StartPeriods must match."))
  }

  PdLabels <- tibble(Labels=OrderedPds, Pd.Num=1:J)
  Start_js <- tibble(Clusters=Clusters,
                     Labels=StartPeriods)
  if (sum(StartPeriods %in% OrderedPds) < length(StartPeriods)) {
    if (is.numeric(OrderedPds) & is.numeric(StartPeriods) &
        sum(OrderedPds[order(OrderedPds)]==OrderedPds)==J) {
      warning(simpleWarning("Since OrderedPds is numeric, StartPeriods outside this range are treated numerically. NAs/Infs are considered never-treated."))
      min_Ord <- min(OrderedPds)
      Start_js <- Start_js %>%
        dplyr::mutate(Pd.Num=if_else(is.na(Labels) | is.infinite(Labels),
                                     Inf,
                                     Labels-min_Ord+1))
    } else {
      warning(simpleWarning("All StartPeriods not in OrderedPds are considered never-treated units within this setting."))
      Start_js <- Start_js %>%
        left_join(PdLabels, by=join_by(Labels)) %>%
        dplyr::mutate(Pd.Num=if_else(is.na(Pd.Num),Inf,Pd.Num))
    }
  } else {
    Start_js <- Start_js %>%
      dplyr::left_join(PdLabels, by=join_by(Labels))
  }

  Start_js <- Start_js %>%
    dplyr::arrange(Pd.Num) %>%
    dplyr::mutate(Cl.Num=1:N) %>%
    dplyr::rename(Start_j = Pd.Num) %>%
    dplyr::select(Clusters,Cl.Num,Start_j)

  return(list(Start_js=Start_js,
              PdLabels=PdLabels %>% dplyr::rename(Periods=Labels) %>%
                dplyr::select(Pd.Num,Periods),
              N=N,
              J=J))
}

# gen_A function

## Helper Functions

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

## Main gen_A function
### Inputs: N (number of periods), J (number of units)
### Output: A DF that converts observations y to two-by-two DID estimators d

gen_A <- function(N,J) {
  J2 <- J*(J-1)/2
  Adot <- gen_Adot(J)

  do.call(rbind,
          lapply(X=1:(N-1),
                 FUN=function(n1) do.call(rbind,
                                          lapply(X=(n1+1):N,
                                                 FUN=function(n2) gen_Arow(N,J,n1,n2,J2,Adot)))))
}

# gen_D function

## Main gen_D function

### Inputs: N (number of periods), J (number of units)
### Output: A DF where each row corresponds to an entry in d,
#### giving its i, i', j, j' (unit and period) indices

gen_D <- function(N,J) {
    Cl.df <- data.frame(i=unlist(sapply(1:(N-1), FUN=function(x) rep(x,N-x))),
                        i.prime=unlist(sapply(1:(N-1), FUN=function(x) (x+1):N)))
    Pd.df <- data.frame(j=unlist(sapply(1:(J-1), FUN=function(x) rep(x,J-x))),
                        j.prime=unlist(sapply(1:(J-1), FUN=function(x) (x+1):J)))
    return(dplyr::cross_join(Cl.df, Pd.df))
}

# gen_Theta function

## Main gen_Theta function

### Inputs:
#### js_Obj: Output of gen_js function
#### Assumption: Number 1--5 corresponding to the Assumption Settings
### Output: List of the following:
#### All: DF of the unique treatment effects
#### Full: DF of all treated cells with their corresponding treatment effect
#### Schematic: a matrix of a stylized stepped-wedge diagram showing the theta number for each treated cell

gen_Theta <- function(js_Obj,Assumption) {
  Theta <- js_Obj$Start_js %>%
    dplyr::cross_join(js_Obj$PdLabels) %>%
    dplyr::filter(Pd.Num >= Start_j) %>%
    mutate(Diff=Pd.Num - Start_j + 1)

  if (Assumption==5) {
    JoinBy <- NULL
  } else if (Assumption==4) {
    JoinBy <- c("Pd.Num")
  } else if (Assumption==3) {
    JoinBy <- c("Diff")
  } else if (Assumption==2) {
    JoinBy=c("Pd.Num","Diff")
  } else if (Assumption==1) {
    JoinBy=c("Cl.Num","Pd.Num","Diff")
  } else {
    stop(simpleError("Assumption must be a value 1 through 5 corresponding to the desired assumption setting."))
  }

  if(is.null(JoinBy)) {
    All <- tibble(Theta=1)
  } else {
    All <- Theta %>% dplyr::select(all_of(JoinBy)) %>%
      dplyr::arrange(across(JoinBy)) %>%
      dplyr::distinct() %>%
      dplyr::mutate(Theta=row_number())
  }
  Full <- merge(Theta, All, by=JoinBy)

  Schematic <- matrix(data=0,
                      nrow=js_Obj$N,
                      ncol=js_Obj$J)
  Schematic[as.matrix(Full %>% dplyr::select(Cl.Num,Pd.Num))] <- Full$Theta
  dimnames(Schematic) <- list(js_Obj$Start_js$Clusters,js_Obj$PdLabels$Periods)

  return(list(All=All,
              Full=Full,
              Schematic=Schematic))
}

# gen_D function

## Main gen_D function

### Inputs:
#### D: the output DF from gen_D
#### Theta: the output list from gen_Theta
### Output: List of the following:
#### D_aug: the DF D augmented with the numbers of each relevant theta for each row
#### F_mat: the matrix that, when left-multiplied by Theta, gives the expectation of d

gen_F <- function(D,Theta) {
  ThetaF <- Theta$Full %>% dplyr::select(Cl.Num,Pd.Num,Theta)
  D_aug <- D %>% dplyr::select(i,i.prime,j,j.prime,Type,TypeLabel) %>%
    dplyr::left_join(ThetaF %>%
                       dplyr::rename(Pos.i=Theta),
                     by=join_by(i==Cl.Num,
                                j.prime==Pd.Num)) %>%
    dplyr::left_join(ThetaF %>%
                       dplyr::rename(Neg.i=Theta),
                     by=join_by(i==Cl.Num,
                                j==Pd.Num)) %>%
    dplyr::left_join(ThetaF %>%
                       dplyr::rename(Neg.i.prime=Theta),
                     by=join_by(i.prime==Cl.Num,
                                j.prime==Pd.Num)) %>%
    dplyr::left_join(ThetaF %>%
                       dplyr::rename(Pos.i.prime=Theta),
                     by=join_by(i.prime==Cl.Num,
                                j==Pd.Num)) %>%
    mutate(Row=row_number())

  F_mat <- matrix(0,nrow=dim(D_aug)[1], ncol=length(Theta$All$Theta))
  F_mat[as.matrix(D_aug %>% dplyr::select(Row,Pos.i) %>% dplyr::filter(!is.na(Pos.i)))] <-
    F_mat[as.matrix(D_aug %>% dplyr::select(Row,Pos.i) %>% dplyr::filter(!is.na(Pos.i)))]+1
  F_mat[as.matrix(D_aug %>% dplyr::select(Row,Neg.i) %>% dplyr::filter(!is.na(Neg.i)))] <-
    F_mat[as.matrix(D_aug %>% dplyr::select(Row,Neg.i) %>% dplyr::filter(!is.na(Neg.i)))]-1
  F_mat[as.matrix(D_aug %>% dplyr::select(Row,Neg.i.prime) %>% dplyr::filter(!is.na(Neg.i.prime)))] <-
    F_mat[as.matrix(D_aug %>% dplyr::select(Row,Neg.i.prime) %>% dplyr::filter(!is.na(Neg.i.prime)))]-1
  F_mat[as.matrix(D_aug %>% dplyr::select(Row,Pos.i.prime) %>% dplyr::filter(!is.na(Pos.i.prime)))] <-
    F_mat[as.matrix(D_aug %>% dplyr::select(Row,Pos.i.prime) %>% dplyr::filter(!is.na(Pos.i.prime)))]+1
  return(list(D_aug=D_aug, F_mat=F_mat))
}

# gen_ADFT function

## Main gen_ADFT function

### Inputs:
#### Clusters, StartPeriods, OrderedPds (as in gen_js input)
#### Assumption (as in gen_Theta input)
### Output: List of the following:
#### A_mat: output from gen_A
#### D_aug, F_mat from gen_F
#### Theta from gen_Theta
#### N, J from gen_js

gen_ADFT <- function(Clusters, StartPeriods, OrderedPds,
                    Assumption=5) {
  js_Obj <- gen_js(Clusters, StartPeriods, OrderedPds)

  A_mat <- gen_A(N=js_Obj$N, J=js_Obj$J)

  Theta <- gen_Theta(js_Obj, Assumption)

  D <- gen_D(N=js_Obj$N,J=js_Obj$J) %>%
    left_join(js_Obj$Start_js %>%
                dplyr::rename(i=Cl.Num,Start.i=Start_j,
                              Cl.i=Clusters), by="i") %>%
    left_join(js_Obj$Start_js %>%
                dplyr::rename(i.prime=Cl.Num,Start.i.prime=Start_j,
                              Cl.i.prime=Clusters), by="i.prime") %>%
    mutate(i.j.leadlag=j-Start.i+1,
           ip.j.leadlag=j-Start.i.prime+1,
           i.jp.leadlag=j.prime-Start.i+1,
           ip.jp.leadlag=j.prime-Start.i.prime+1,
           i.type=if_else(i.j.leadlag > 0,"Always-Treated",
                          if_else(i.jp.leadlag <= 0, "Always-Control","Switch")),
           ip.type=if_else(ip.j.leadlag > 0,"Always-Treated",
                           if_else(ip.jp.leadlag <= 0, "Always-Control","Switch")),
           Type=if_else(i.type=="Always-Control",1,
                        if_else(i.type=="Switch",
                                if_else(ip.type=="Switch",4,2),
                                if_else(ip.type=="Always-Control",3,
                                        if_else(ip.type=="Switch",5,6)))),
           TypeLabel=factor(Type,
                            levels=1:6, labels=c("Both Always-Control",
                                                 "Switch vs. Always-Control",
                                                 "Always-Treated vs. Always-Control",
                                                 "Both Switch",
                                                 "Always-Treated vs. Switch",
                                                 "Both Always-Treated")))

  Fres <- gen_F(D, Theta)
  return(list(A_mat=A_mat, D_aug=Fres$D_aug, F_mat=Fres$F_mat, Theta=Theta,
              N=js_Obj$N, J=js_Obj$J))
}

# create_V function

## main create_V function

### Inputs:
#### Length: the number of unique Theta values
#### NonZero: the indices that should be non-zero
#### Values: the weights for non-zero indices, in the same order as NonZero
####  (if Null, defaults to equal weight, summing to 1)
### Output:
#### A vector of length Length with the corresponding values, to be used as estimand weights

create_V <- function(
    Length,
    NonZero,
    Values = NULL) {
  v <- rep(0, Length)
  if (is.null(Values)) {
    v[NonZero] <- 1 / length(NonZero)
  } else {
    v[NonZero] <- Values
  }
  return(v)
}

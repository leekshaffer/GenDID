#######################################
###### File: Var_Min.R ################
#######################################

# Helper functions

## opt_fn function

### Inputs:
#### x: the vector, which with a preceding 1, is on the outsides of the quadratic form
#### Inner: the square matrix which is inside the quadratic form, of dimension 1 more than x
### Output: matrix resulting from quadratic form multiplication

opt_fn <- function(x,Inner) {
  matrix(data=c(1,x), nrow=1) %*% Inner %*% matrix(data=c(1,x), ncol=1)
}

## min_var_single function

### Inputs:
#### base.w: a single base observations weight that is a solution to the constraint
#### Add.Obs.w: the additional observation weight bases that do not affect the expectation
#### ASigAT: the matrix resulting from the quadratic form of A and Sigma
#### method: the method for optimization in the optim() function (default="CG")
#### maxit: passed to the control input of optim() (default=10000)
### Outputs: List of the following:
#### Cvec: the optimal observation weights
#### Variance: the optimal (minimum) variance, based on input Sigma
#### DID.weights: one set of corresponding DID weights

min_var_single <- function(base.w,Add.Obs.w,ASigAT,method="CG", maxit=10000) {
  W_mat <- as.matrix(cbind(as.matrix(base.w, ncol=1),Add.Obs.w))
  colnames(W_mat) <- c("base.w",colnames(Add.Obs.w))
  Inn_mat <- t(W_mat) %*% ASigAT %*% W_mat
  if (ncol(Add.Obs.w)==1) { ## Solved quadratic case for k=1:
    if (Inn_mat[2,2] > 0) {
      Min <- -1*Inn_mat[1,2]/Inn_mat[2,2]
      Var <- Inn_mat[1,1]+2*Min*Inn_mat[1,2]+Min^2*Inn_mat[2,2]
      # Var <- as.numeric(opt_fn(x=Min, Inner=Inn_mat))
    } else {
      stop(simpleError("No Minimum Found."))
    }
  } else {
    Opt <- optim(par=rep(1,ncol(Add.Obs.w)),
                 fn=opt_fn,
                 method=method, control=list(maxit=maxit),
                 Inner=Inn_mat)
    Min <- Opt$par
    Var <- Opt$value
  }
  return(list(Cvec=c(1,Min),
              Variance=Var,
              DID.weights=W_mat %*% matrix(c(1,Min),ncol=1)))
}

# min_var function

### Inputs:
#### SolveOut: The output from the Solve_Assumption function
#### Sigma: The variance matrix to use for variance minimization
#### Drop_ADFT: If TRUE, the ADFT will not be passed through to the output
#### method: the method for optimization in the optim() function (default="CG")
#### maxit: passed to the control input of optim() (default=10000)
### Outputs: List of the following:
#### Cvec: the optimal observation weights
#### Variance: the optimal (minimum) variance
#### DID.weights: one set of corresponding DID weights
#### ADFT: the ADFT object passed through from SolveOut

## A wrapper that outputs the full results of the minimization and can handle a matrix of v's:
## Note the Var that is returned is based on the input Sigma,
### and should generally be used only for relative comparisons of weightings
min_var <- function(SolveOut,
                    Sigma,
                    Drop_ADFT=FALSE,
                    method="CG",maxit=10000) {
  A_mat <- SolveOut$ADFT$A_mat
  solve_obj <- SolveOut$Solve
  if (nrow(Sigma) != ncol(Sigma) | nrow(Sigma) != ncol(A_mat)) {
    stop(simpleError("Sigma must be a square matrix with the same number of columns as A."))
  }

  base.w <- solve_obj$DID.weights %>% dplyr::select(starts_with("w.base")) %>%
    dplyr::rename_all(~stringr::str_replace(.,"^w.base.",""))
  Add.Obs.w <- solve_obj$DID.weights %>% dplyr::select(starts_with("Add.Obs.weights"))
  ASigAT <- A_mat %*% Sigma %*% t(A_mat)

  if (ncol(Add.Obs.w) == 0) {
    print("There is a unique solution to the constraint.")
    res <- list(Cvec=c(1),
                Variance=diag(t(base.w) %*% ASigAT %*% as.matrix(base.w)),
                DID.weights=base.w,
                Obs.weights=solve_obj$Obs.weights %>% dplyr::select(starts_with("ATw.base")))
  } else {
    if (ncol(base.w)==1) { ## A single target estimand
      res1 <- min_var_single(base.w,Add.Obs.w,ASigAT,
                            method=method,maxit=maxit)
      Cvec <- as.matrix(res1$Cvec, ncol=1)
      colnames(Cvec) <- colnames(base.w)
      Var <- as.matrix(res1$Variance, ncol=1)
      colnames(Var) <- colnames(base.w)
      DID.weights <- as.matrix(res1$DID.weights, ncol=1)
      colnames(DID.weights) <- colnames(base.w)
      res <- list(Cvec=Cvec,
                  Variance=Var,
                  DID.weights=DID.weights,
                  Obs.weights=t(A_mat) %*% DID.weights)
    } else { ## Multiple target estimands
      res1 <- apply(base.w, MARGIN=2,
                   FUN=function(col) min_var_single(col,Add.Obs.w,ASigAT,
                                                    method=method,maxit=maxit))
      DID.weights <- do.call(cbind, lapply(res1, FUN=function(x) x$DID.weights))
      colnames(DID.weights) <- colnames(base.w)
      res <- list(Cvec=do.call(cbind, lapply(res1, FUN=function(x) x$Cvec)),
                  Variance=do.call(cbind, lapply(res1, FUN=function(x) x$Variance)),
                  DID.weights=DID.weights,
                  Obs.weights=t(A_mat) %*% DID.weights)
    }
  }
  if (Drop_ADFT) {
    return(res)
  } else {
    return(list(ADFT=SolveOut$ADFT,
                MV=res))
  }
}

# Methods for nested comparison estimates and weights:

# Comp_Ests function

## main Comp_Ests function

### Inputs:
#### ADFT_obj: output from gen_ADFT
#### estimator: A character vector of estimators to compute. Options are:
####  "TW": Two-way fixed effects estimator; see Goodman-Bacon (2021).
####  "CS": Callaway and Sant'Anna estimator.
####  "SA": Sun and Abraham estimator.
####  "CH": Chaisemartin and D’Haultfœuille estimator.
####  "CO": Crossover-based estimator; see Kennedy-Shaffer et al. (2020).
### Output: List of the following:
#### Weights: A matrix containing the computed weights for the specified estimators.
#### Full: A detailed list containing intermediate and final results for each estimator.

Comp_Ests <- function(ADFT_obj,
                      estimator=c("TW","CS","SA","CH","CO")) {
  ### First, augment D with start times:
  Starts <- ADFT_obj$Theta$Full %>%
    dplyr::select(Cl.Num,Start_j) %>%
    distinct()
  if (max(Starts$Cl.Num) < ADFT_obj$N) {
    Starts <- Starts %>%
      bind_rows(tibble(Cl.Num=(max(Starts$Cl.Num)+1):ADFT_obj$N,
                       Start_j=Inf))
  }
  D_use <- ADFT_obj$D_aug %>% dplyr::select(i,i.prime,j,j.prime,Type,TypeLabel) %>%
    left_join(Starts %>% rename(i.start=Start_j),
              by=join_by(i == Cl.Num)) %>%
    left_join(Starts %>% rename(i.prime.start=Start_j),
              by=join_by(i.prime == Cl.Num))

  res <- NULL

  if ("TW" %in% estimator) { ## Based on weights from Goodman-Bacon (2021):
    Starts_n <- Starts %>% group_by(Start_j) %>%
      dplyr::summarize(n.Clust=n()) %>%
      dplyr::mutate(n=n.Clust/sum(n.Clust),
                    Dbar=if_else(is.infinite(Start_j),0,
                                 (ADFT_obj$J-Start_j+1)/ADFT_obj$J)) %>%
      dplyr::rename(G=Start_j)
    Cross_starts <- cross_join(Starts_n, Starts_n,
                               suffix=c(".k",".l")) %>%
      dplyr::filter(G.k != G.l, Dbar.k > 0)
    TW.Weights <- Cross_starts %>%
      dplyr::mutate(s.num=if_else(Dbar.l==0, n.k*(n.l)*Dbar.k*(1-Dbar.k),
                                  if_else(G.k < G.l, n.k*n.l*(Dbar.k-Dbar.l)*(1-Dbar.k),
                                          -1*n.k*n.l*Dbar.k*(Dbar.l-Dbar.k))),
                    n.DID=if_else(Dbar.l==0, (ADFT_obj$J-G.k+1)*(G.k-1)*n.Clust.k*n.Clust.l,
                                  if_else(G.k < G.l, (G.l-G.k)*(G.k-1)*n.Clust.k*n.Clust.l,
                                          (ADFT_obj$J-G.k+1)*(G.k-G.l)*n.Clust.k*n.Clust.l))) %>%
      dplyr::mutate(s=s.num/sum(abs(s.num)),
                    W_TW=s/n.DID)
    TW_w <- D_use %>% mutate(G.k=if_else(Type==2, i.start,
                                         if_else(Type==5, i.prime.start, NA)),
                             G.l=if_else(Type==2, i.prime.start,
                                         if_else(Type==5, i.start, NA))) %>%
      left_join(TW.Weights %>% dplyr::select(G.k,G.l,W_TW),
                by=c("G.k","G.l")) %>%
      dplyr::mutate(W_TW=if_else(is.na(W_TW), 0, W_TW))
    res <- c(res, list(W_TW=TW_w %>% dplyr::select(W_TW)))
  }

  if ("CS" %in% estimator) {
    Group <- unique(Starts$Start_j[Starts$Start_j < max(Starts$Start_j)])
    Time <- unique(D_use %>% filter(Type %in% 1:3) %>% pull(j.prime))
    CS_key <- crossing(Group,Time) %>% mutate(Col=1:n())
    CS_gt_w_int <- apply(as.matrix(CS_key),
                         MARGIN=1,
                         FUN=function(row) CS_wt_fun(D_use, row["Group"], row["Time"]))
    CS_gt_w <- apply(CS_gt_w_int, 2, FUN=function(x) x/sum(abs(x)))
    CS_key_agg <- CS_key %>% mutate(Type=if_else(Time >= Group,"Post","Pre"),
                                    EventTime=Time-Group) %>%
      left_join(Starts %>% group_by(Start_j) %>% dplyr::summarize(Num=n()),
                by=join_by(Group==Start_j))
    CS_key_agg <- CS_key_agg %>%
      left_join(CS_key_agg %>% group_by(EventTime) %>% dplyr::summarize(Num_dynamic=sum(Num))) %>%
      left_join(CS_key_agg %>% group_by(Group) %>%
                  dplyr::summarize(Num_group=sum(if_else(Type=="Post",Num,0)))) %>%
      left_join(CS_key_agg %>% group_by(Time) %>%
                  dplyr::summarize(Num_calendar=sum(if_else(Type=="Post",Num,0)))) %>%
      mutate(W_simple=if_else(Type=="Post",Num,0)/sum(if_else(Type=="Post",Num,0)),
             W_dynamic=(if_else(Type=="Post",Num/Num_dynamic,0)/sum(if_else(Type=="Post",Num/Num_dynamic,0))),
             W_group=(if_else(Type=="Post",Num/Num_group,0)/sum(if_else(Type=="Post",Num/Num_group,0))),
             W_calendar=(if_else(Type=="Post",Num/Num_calendar,0)/sum(if_else(Type=="Post",Num/Num_calendar,0))))
    CS_w <- CS_gt_w %*% as.matrix(CS_key_agg %>% dplyr::select(starts_with("W_")))
    res <- c(res, list(Key_CS=CS_key_agg, W_CS_gt=CS_gt_w, W_CS=CS_w))
  }

  if ("SA" %in% estimator) {
    Rel_Pds <- unique(D_use %>% filter(Type %in% 1:3) %>%
                        mutate(RelPd=j.prime-i.start) %>% pull(RelPd))
    Rel_Pds <- sort(Rel_Pds[Rel_Pds != -1])
    Ctrl.is <- Starts %>% filter(Start_j==max(Start_j)) %>% pull(Cl.Num)
    CATT_el_key <- crossing(Period=Rel_Pds,
                            Cohort=Starts %>% filter(Start_j != max(Start_j)) %>%
                              pull(Start_j))
    CATT_el_key$Num <- apply(as.matrix(CATT_el_key), MARGIN=1,
                             FUN=function(row) nrow(D_use %>% filter(Type %in% 1:3,
                                                                     i.start==row["Cohort"],
                                                                     i.prime %in% Ctrl.is,
                                                                     j.prime-i.start==row["Period"],
                                                                     j-i.start==-1)) +
                               nrow(D_use %>% filter(Type %in% 1:3,
                                                     i.start==row["Cohort"],
                                                     i.prime %in% Ctrl.is,
                                                     j.prime-i.start==-1,
                                                     j-i.start==row["Period"])))
    CATT_el_key <- CATT_el_key %>% filter(Num != 0) %>%
      mutate(Column=1:n())
    CATT_el_w_int <- apply(CATT_el_key, MARGIN=1,
                           FUN=function(row) SA_wt_fun(D_use, Ctrl.is,
                                                       row["Period"],row["Cohort"]))
    CATT_el_w <- apply(CATT_el_w_int, MARGIN=2,
                       FUN=function(x) x/sum(abs(x)))
    SA_l_key <- CATT_el_key %>% dplyr::select(Period,Num) %>% group_by(Period) %>%
      dplyr::summarize(Num=sum(Num)) %>% mutate(Column=1:n()) %>%
      mutate(Type=if_else(Period >= 0,"Post","Pre"),
             W_SA.W_ATT=if_else(Type=="Post",Num,0)/sum(if_else(Type=="Post",Num,0)))
    SA_l_w <- apply(SA_l_key,
                    MARGIN=1,
                    FUN=function(row) SA_CATT_agg(CATT_el_key,
                                                  CATT_el_w,
                                                  as.numeric(row["Period"]))/as.numeric(row["Num"]))
    res <- c(res,
             list(Key_SA_l=SA_l_key, W_SA_l=SA_l_w,
                  W_SA=SA_l_w %*% as.matrix(SA_l_key %>% dplyr::select(W_SA.W_ATT))))
  }

  if ("CH" %in% estimator) {
    CH_Pt_key <- D_use %>% dplyr::filter(Type==2,
                                         j.prime-j==1,
                                         j.prime==i.start) %>%
      group_by(j.prime) %>%
      dplyr::summarize(Num=n(),
                       Num.Trt=length(unique(i))) %>%
      rename(Period=j.prime) %>%
      mutate(Column=1:n())
    CH_Pt_w_int <- apply(CH_Pt_key, MARGIN=1,
                         FUN=function(row) CH_wt_fun(D_use,row["Period"]))
    CH_Pt_w <- apply(CH_Pt_w_int, 2, function(x) x/sum(abs(x)))
    res <- c(res, list(W_CH=tibble(W_CH.W_M = CH_Pt_w %*% CH_Pt_key$Num.Trt/sum(CH_Pt_key$Num.Trt))))
  }

  if ("CO" %in% estimator) {
    CO_Pt_key <- D_use %>% dplyr::filter(Type %in% c(2,5),
                                         j.prime-j==1,
                                         j.prime==i.start | j.prime==i.prime.start) %>%
      group_by(j.prime,Type) %>%
      dplyr::summarize(Num.Comps=n(),
                       Num.i=length(unique(i)),
                       Num.i.prime=length(unique(i.prime))) %>%
      pivot_wider(id_cols="j.prime",
                  names_from="Type",
                  values_from=starts_with("Num"),
                  names_sep=".",
                  values_fill=0) %>%
      mutate(Num.Comps=Num.Comps.2+Num.Comps.5,
             Num.Trt=max(Num.i.2,Num.i.prime.5),
             Num.Ctrl=Num.i.prime.2+Num.i.5,
             Num.Comps.Res=Num.Comps.2,
             Num.Trt.Res=Num.i.2,
             Num.Ctrl.Res=Num.i.prime.2) %>%
      ungroup() %>%
      dplyr::select(j.prime,Num.Comps,Num.Trt,Num.Ctrl,
                    Num.Comps.Res,Num.Trt.Res,Num.Ctrl.Res) %>%
      rename(Period=j.prime) %>%
      mutate(Column=1:n())
    CO_Pt_w_int <- apply(CO_Pt_key, MARGIN=1,
                         FUN=function(row) CO_wt_fun(D_use,row["Period"]))
    CO_Pt_All <- apply(CO_Pt_w_int, 2,
                       function(x) x/if_else(sum(abs(x))==0,1,sum(abs(x))))
    CO_Pt_All <- apply(CO_Pt_w_int,
                       2, function(x) x/sum(abs(x)))
    CO_Pt_Res <- apply(CO_Pt_w_int,
                       2, function(x) if_else(x<0,0,x)/if_else(sum(x[x>=0])==0,1,sum(x[x>=0])))
    CO_Pt_key <- CO_Pt_key %>%
      mutate(W_CO1=if_else(Num.Comps.Res > 0,1,0)/sum(if_else(Num.Comps.Res > 0,1,0)),
             W_CO2=if_else(Num.Comps.Res > 0,
                           2/(1/Num.Trt.Res+1/Num.Ctrl.Res),
                           0)/sum(if_else(Num.Comps.Res > 0,
                                          2/(1/Num.Trt.Res+1/Num.Ctrl.Res),
                                          0)),
             W_CO3=if_else(Num.Comps > 0,1,0)/sum(if_else(Num.Comps > 0,1,0)))
    CO_w <- bind_cols(CO_Pt_Res %*% as.matrix(CO_Pt_key %>% dplyr::select(c("W_CO1","W_CO2"))),
                      CO_Pt_All %*% as.matrix(CO_Pt_key %>% dplyr::select("W_CO3")))
    res <- c(res,list(Key_CO=CO_Pt_key,
                      W_CO_Pt_All=CO_Pt_All,
                      W_CO_Pt_Res=CO_Pt_Res,
                      W_CO=CO_w))
  }
  Weights <- do.call(cbind, res[paste0("W_",estimator)])
  return(list(Weights=Weights,Full=res))
}

# Comp_Ests_Weights function

## main Comp_Ests_Weights function

### Inputs:
#### ADFT_obj: output from gen_ADFT
#### estimator: A character vector of estimators to compute; see Comp_Ests for options
### Output: List of the following:
#### DID.weights:  A matrix containing computed DID weights for the specified estimators.
#### Obs.weights: A matrix containing computed observation weights for the specified estimators.

Comp_Ests_Weights <- function(ADFT_obj,
                              estimator=c("TW","CS","SA","CH","CO","NP")) {
  DID.weights <- NULL
  Obs.weights <- NULL
  if (sum(c("TW","CS","SA","CH","CO") %in% estimator) > 0) {
    DID.ests <- estimator[estimator %in% c("TW","CS","SA","CH","CO")]
    DID.weights <- bind_cols(DID.weights,Comp_Ests(ADFT_obj,
                                                   DID.ests)$Weights)
    Obs.weights <- bind_cols(Obs.weights, t(ADFT_obj$A_mat) %*% as.matrix(DID.weights))
  }
  if ("NP" %in% estimator) {
    N_Trt <- apply(ADFT_obj$Theta$Schematic, 2, function(x) sum(x > 0))
    N_Ctrl <- ADFT_obj$N - N_Trt
    W_Trt <- if_else(N_Trt==0 | N_Ctrl==0,0,1/N_Trt)
    W_Ctrl <- if_else(N_Trt==0 | N_Ctrl==0,0,-1/(N_Ctrl))
    Mat <- matrix(c(W_Trt, N_Trt, W_Ctrl, N_Ctrl), nrow=4, byrow=TRUE)
    W_Mat <- apply(Mat, 2, FUN=function(col) c(rep(col[1], col[2]), rep(col[3], col[4])))
    NP_w_int <- as.vector(t(W_Mat))
    Mult_Eq <- if_else(N_Trt > 0 & N_Ctrl > 0, 1, 0)/sum(if_else(N_Trt > 0 & N_Ctrl > 0, 1, 0))
    Mult_ATT <- if_else(N_Trt > 0 & N_Ctrl > 0, N_Trt, 0)/sum(if_else(N_Trt > 0 & N_Ctrl > 0, N_Trt, 0))
    Mult_IV <- if_else(N_Trt > 0 & N_Ctrl > 0, (1/N_Trt + 1/N_Ctrl)^(-1),0)/sum(if_else(N_Trt > 0 & N_Ctrl > 0, (1/N_Trt + 1/N_Ctrl)^(-1),0))
    Obs.weights <- bind_cols(Obs.weights,
                             W_NP_Eq=NP_w_int * rep(Mult_Eq, times=ADFT_obj$N),
                             W_NP_ATT=NP_w_int * rep(Mult_ATT, times=ADFT_obj$N),
                             W_NP_IV=NP_w_int * rep(Mult_IV, times=ADFT_obj$N))
  }
  return(list(DID.weights=DID.weights,
              Obs.weights=Obs.weights))
}

# Helper functions for specific methods

## CS_wt_fun function:

CS_wt_fun <- function(D_use,group,time) {
  if (time >= group) {
    if_else(D_use$Type %in% 1:3 & D_use$i.start==group & D_use$i.prime.start != group &
              D_use$j.prime==time & D_use$j==group-1,
            1,0)
  } else {
    if_else(D_use$Type & D_use$j.prime & D_use$j==time-1,
            if_else(D_use$i.start==group & D_use$i.prime.start != group, 1,
                    if_else(D_use$i.start != group & D_use$i.prime.start == group, -1, 0)), 0)
  }
}

## SA_wt_fun function:

SA_wt_fun <- function(D_use,Ctrl.is,period,cohort) {
  if (period >= 0) {
    if_else(D_use$Type %in% 1:3 & D_use$i.start==cohort &
              D_use$i.prime %in% Ctrl.is &
              as.numeric(D_use$j.prime) - as.numeric(D_use$i.start) == period &
              as.numeric(D_use$j) - as.numeric(D_use$i.start) == -1,
            1, 0)
  } else {
    if_else(D_use$Type %in% 1:3 & D_use$i.start==cohort &
              D_use$i.prime %in% Ctrl.is &
              as.numeric(D_use$j.prime) - as.numeric(D_use$i.start) == -1 &
              as.numeric(D_use$j) - as.numeric(D_use$i.start) == period,
            -1, 0)
  }
}

## SA_CATT_agg function:

SA_CATT_agg <- function(key_tbl, wt_tbl, period) {
  key_vals <- key_tbl %>% dplyr::filter(Period==period)
  wt_tbl[,key_vals$Column,drop=FALSE] %*% key_vals$Num
}

## CH_wt_fun function:

CH_wt_fun <- function(D_use,period) {
  if_else(D_use$Type==2 & D_use$j.prime==D_use$i.start &
            D_use$j.prime==period & D_use$j == period-1,
          1,0)
}

## CO_wt_fun function:

CO_wt_fun <- function(D_use,period) {
  if_else(D_use$Type==2,
          if_else(D_use$j.prime==D_use$i.start &
                    D_use$j.prime==period &
                    D_use$j==period-1, 1, 0),
          if_else(D_use$Type==5,
                  if_else(D_use$j.prime == D_use$i.prime.start &
                            D_use$j.prime == period & D_use$j == period-1,
                          -1,0),0))
}

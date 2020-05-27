#include <RcppArmadillo.h>
#include <vector>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]


/*Calculate the grouplasso portion of the loss function*/
// [[Rcpp::export]]
double grouplassoLossCpp(const arma::uvec& groups, const arma::vec& grouplasso, const arma::mat& W, const int& Q){
    if(arma::accu(grouplasso) == 0){
        return 0;
    } else {
        
        double out = 0;
        for(int j = 0; j < Q; j++){

            unsigned int from = 0;
            int to = -1;

            for(arma::uword i = 0; i < groups.n_elem; i++){ 
                to += groups[i];
                out += std::sqrt(groups[i]) * arma::norm(W.col(j).subvec(from, to)) * grouplasso[j];
                from += groups[i];
                
                } 
            }
        return out; 
    }
}



/*Calculate the grouplasso portion of the loss function
  note this is a void function*/
// [[Rcpp::export]]
void grouplassoPenaltyCpp(const arma::uvec& groups, const arma::mat& W, arma::mat& DG, const int& Q){
    for(int j = 0; j < Q; j++){

        unsigned int from = 0;
        int to = -1;
 
        for(arma::uword i = 0; i < groups.n_elem; i++){ 

            to += groups[i];
            double filler = std::sqrt(groups(i))/2 * 1/arma::norm(W.col(j).subvec(from, to));

            DG.col(j).subvec(from, to).fill(filler); 
            
            from += groups[i];
        }
    }
    //return DG;
}


//calculate the elitist part of the loss function
// [[Rcpp::export]]
double elitistLossCpp(const arma::uvec& groups, const arma::vec& elitistlasso, const arma::mat& W, const int& Q){
    if(arma::accu(elitistlasso) == 0){
        return 0;
    } else {
        double out = 0;
        for(int j = 0; j < Q; j++){

            unsigned int from = 0;
            int to = -1;

            for(arma::uword i = 0; i < groups.n_elem; i++){ 

                to += groups[i];
                out += std::pow(arma::norm(W.col(j).subvec(from, to), 1), 2) * elitistlasso[j];
                from += groups[i];

            } 
        }
        return out; 
    }
}



/*Calculate the elist part of the penalty, note this function is a void function,
 also note that if all elements in W are constraints to zero, DE/arma::abs(W) will be 0/0 == NaN,
 in the coordinate descent step NaN gets multiplied with 0 leading to an error, therefore if DE contains exact zeroes
 which I assume can only be the case if the user supplied them, I will change the values into ones, to prevent an error
 in the coordinate descent step */
// [[Rcpp::export]]
void elitistPenaltyCpp(const arma::uvec& groups, const arma::mat& W, arma::mat& DE, const int& Q){
    for(int j = 0; j < Q; j++){

        unsigned int from = 0;
        int to = -1;
 
        for(arma::uword i = 0; i < groups.n_elem; i++){ 

            to += groups[i];
            double filler = arma::norm(W.col(j).subvec(from, to), 1);

            DE.col(j).subvec(from, to).fill(filler); 
            
            from += groups[i];
        }
    }

    DE.elem(find(DE == 0)).ones();
    DE =  DE / arma::abs(W);
}

//Calculates the loss function
// [[Rcpp::export]]
long double lossFunctionCpp(const arma::mat& X, const arma::mat& W, const arma::mat& P, const arma::uvec& groups, const arma::vec& ridge, const arma::vec& lasso, const arma::vec& grouplasso, const arma::vec& elitistlasso, const double& Q){

        double out = arma::accu(arma::pow((X - X * W * P.t()), 2))
        +  elitistLossCpp(groups, elitistlasso, W, Q)
        +  grouplassoLossCpp(groups, grouplasso, W, Q)
        +  arma::accu(arma::pow(W, 2)*arma::diagmat(ridge))
        +  arma::accu(arma::abs(W)*arma::diagmat(lasso));

        return out;
} 


//Updates the elements of W once component after componet, could be done in parallel
// [[Rcpp::export]]
void coordinateDescentStepCpp(const arma::mat& XtXP, const arma::mat& XtX, arma::mat& W, const arma::mat& constraints, const int& Q, const int& J, const arma::mat& D){
    double num, denom;
    for(int q = 0; q < Q; q++){
        for(int j = 0; j < J; j++){
            if(constraints(j, q) != 0){
                if(j == 0){
                    num = XtXP(j, q) - arma::as_scalar(XtX.row(j).subvec(1, J-1) * W.col(q).subvec(1, J-1));
                } else if (j == J-1){
                    num = arma::as_scalar(XtXP(j, q) - XtX.row(j).subvec(0, J-2) * W.col(q).subvec(0, J-2));
                } else {
                    num = arma::as_scalar(XtXP(j, q) - XtX.row(j).subvec(0, j-1) * W.col(q).subvec(0, j-1)
                        - XtX.row(j).subvec(j+1, J-1) * W.col(q).subvec(j+1, J-1)); 
                }
                denom = XtX(j, j) + D(j, q);
                W(j, q) = num / denom;
            }
        }
    }
}


//Find the mimimum of the majorizing function given P, using coordinate descent 
// [[Rcpp::export]]
void coordinateDescentCpp(const arma::mat& X, const arma::mat& XtX, const arma::mat& XtXP, const int& itrCoor, arma::mat& W, const arma::mat& constraints, const arma::mat& P, const arma::mat& D, const int& Q, const int& J, const bool& printLoss){

    double tol = std::pow(10.0, -8);
    double a, b;
    arma::vec loss = arma::vec(itrCoor + 1); 
    loss(0) = arma::datum::inf;
    for(int i = 0; i < itrCoor; i++){
        a = 0;
        for(int q = 0; q < Q; q++){
            a += arma::as_scalar(arma::pow(W.col(q), 2).t() * D.col(q));
        }

        b = arma::accu(arma::pow(X - X * W * P.t(), 2));
        loss(i + 1) = a + b;
        if(printLoss && (i+1) % 1000 == 0 ){
            Rcpp::Rcout << "iterations coordinate descent step: "  << i+1 << " at loss value: " <<
               loss(i+1) << "\n";
        }
        if(loss(i) - loss(i + 1) < tol){
            return;
        }
        coordinateDescentStepCpp(XtXP, XtX, W, constraints, Q, J, D);
    }
}

//store in a std::vector, uvec's containing the row indices of the free coefficients in W
//note: with no argument supplied, find will return the indices containing the non-zero elements
// [[Rcpp::export]]
std::vector<arma::uvec>  makeConstraints(const arma::mat& constraints, const int& Q){
    
    std::vector<arma::uvec> out;
    out.reserve(Q);

    for(int q = 0; q < Q; q++){
        arma::uvec indices = arma::find(constraints.col(q));
        out.push_back(indices);

    }
    return out;

}

//Full algorithm
// [[Rcpp::export]]
Rcpp::List newAlgoCpp(const arma::mat& X, const arma::vec& ridge, const arma::vec& lasso, const arma::mat& constraints, const arma::vec& grouplasso, const arma::vec& elitistlasso, arma::uvec& groups, const int& Q, const int& itr, arma::mat Wstart, int nStarts = 1, bool printLoss = true, bool coorDec = false) {

    //Object intialization    
    double tol = std::pow(10.0, -6);
    double minLoss, minLossGlobal = arma::datum::inf;
    bool converged = false;
    int J = X.n_cols;
    arma::mat P = arma::randn(J, Q);
    arma::mat DG = arma::randn(J, Q);
    //arma::mat DE = arma::randn(J, Q);
    arma::mat DE = arma::zeros(J, Q);
    arma::mat DR = arma::ones(J, Q);  

    arma::mat U, U2, V, V2, DL, Dsup;
    arma::vec D, D2;
    arma::svd(U, D, V, X);
    arma::mat W = V.cols(0, Q - 1);
    arma::mat XtX = X.t() * X;
    Rcpp::List ret;

    /* needed for the non-contigues submatrix view method: submat*/ 
    std::vector<arma::uvec> indrow = makeConstraints(constraints, Q);
    std::vector<arma::uvec> indcol;
    indcol.reserve(Q);
    
    for(arma::uword q = 0; q < arma::uword(Q); q++){
        arma::uvec a = {q};
        indcol.push_back(a);
    } 

    for(int j = 0; j < nStarts; j++){

        arma::vec loss = arma::vec(itr+1); 
        loss.fill(arma::datum::inf);

        
        /* Rcpp cannot handle matrix intialization therefore
         * the user has to give starting values for the algorithm.
         * if the user supplied a matrix that sums to zero normal
         * initialization follows.         
         * If the user does not want to use custom starting values do:
         */

        if(arma::accu(Wstart) == 0){
            /*the first start of the algorithm will be "warm", the consecutive 
              starts will be random starts */
            //Rcpp::Rcout << "check" << "\n";
            arma::svd(U, D, V, X);
            W = V.cols(0, Q - 1);
            W += arma::randu<arma::mat>(size(W)) * j;
            W.elem(find(constraints == 0)).zeros();
        } else {
            /* If the starting value matrix does not sum to zero,
             * meaning the user wants custom starting values, do:
             */
            W = Wstart;
            W += arma::randu<arma::mat>(size(W)) * j;
            W.elem(find(constraints == 0)).zeros();
        }
        /* start algorithm */
        for(int i = 0; i < itr; i++){
            Rcpp::checkUserInterrupt();

            loss(i + 1) = lossFunctionCpp(X, W, P, groups, ridge, lasso, grouplasso, elitistlasso, Q);

            if(printLoss){
                Rcpp::Rcout << loss(i+1) << "\n";
            }

            //procruste rotation least squares P given W
            arma::svd(U2, D2, V2, XtX * W);
            P = U2.cols(0, Q-1) * V2.t();
            
            //lasso
            DL = arma::pow(arma::abs(W), -1);
            DL.elem(find(DL > std::pow(10.0, 6))).fill(std::pow(10.0, 6));
            //group lasso
            grouplassoPenaltyCpp(groups, W, DG, Q); 
            DG.elem(find(DG > std::pow(10.0, 6))).fill(std::pow(10.0, 6));

            //elitist lasso
            if(arma::accu(elitistlasso) != 0){
                elitistPenaltyCpp(groups, W, DE, Q); 
                DE.elem(find(DE > std::pow(10.0, 6))).fill(std::pow(10.0, 6));
            }

            //scale each column of the diagonal matrices by multiplying with penalties in diagonal matrices 
            Dsup = DL*arma::diagmat(lasso/2) + DR*arma::diagmat(ridge) + DG*arma::diagmat(grouplasso/2)
                + DE*arma::diagmat(elitistlasso);

            arma::mat XtXP = XtX * P;

            /* Given P find solution for W with coordinate descent or
            through inverses*/
            if(coorDec){
                coordinateDescentCpp(X, XtX, XtXP, 1000, W, constraints, P, Dsup, Q, J, printLoss);
            } else {
                for(int q = 0; q < Q; q++){
                    W.submat(indrow[q], indcol[q]) =
                        (arma::diagmat(Dsup.submat(indrow[q], indcol[q])) +
                         XtX.submat(indrow[q], indrow[q])).i() * XtXP.submat(indrow[q], indcol[q]);
                }
            }

            //if converged break
            if(loss(i) - loss(i+1) < tol){
                Rcpp::Rcout << "converged" << "\n";
                converged = true;
                break; 
            }
        }

        loss = loss(find(loss != arma::datum::inf));
        minLoss = loss(loss.size() - 1);       

        if(minLoss < minLossGlobal){

            minLossGlobal = minLoss;
            //set small elements to zero
            W.elem(find(arma::abs(W) < std::pow(10.0, -4))).zeros(); 
            ret["W"] = W;
            ret["P"] = P; 
            ret["loss"] = minLoss;
            ret["converged"] = converged;
        } 
    }

    return ret;
}



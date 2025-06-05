#include "SNPvector.h"
#include "SNPmatrix.h"
#include "LD_EM.h"
#include "richArray.h"
#include <stdexcept>
#include "debug.h"


#ifndef _LD_matrix_
#define _LD_matrix_

enum LDalgorithm { moments, EM };

// NOTE this class is created as a workaround for the impossibility 
// of partial specialization of templated functions
//
// calcule le LD d'une pair de SNP...
template<LDalgorithm Algo, typename scalar_t = double>
class LD_pair_f;

// ... par la méthode des moments
template<typename scalar_t> 
class LD_pair_f<LDalgorithm::moments, scalar_t>{
  public:
  inline scalar_t operator()(SNPmatrix & M, size_t i1, size_t i2) {
    SNPvector & snp1 = *(M.getSNP(i1));
    SNPvector & snp2 = *(M.getSNP(i2));
    return snp1.LD<scalar_t>(snp2);
  }
};

// ... par l'algo EM
template<typename scalar_t>
class LD_pair_f<LDalgorithm::EM, scalar_t> {
  public:
  inline scalar_t operator()(SNPmatrix & M, size_t i1, size_t i2) {
    richArray<9, unsigned int> table;
    SNPvector & snp1 = *(M.getSNP(i1));
    SNPvector & snp2 = *(M.getSNP(i2));
    snp1.contingency(snp2, table);
    return LD_EM<scalar_t>(table);
  }
};




// remplit une matrice carrée de LD des SNPS i avec c1 <= i <= c2 
// doit être ok avec n'importe quelle classe de matrice qui a des membres nrow() ncol() et l'affectation par M(i,j) = ...
// (column major mode matrix)
 
template<LDalgorithm Algo, typename scalar_t = double, typename matrixType>
void LD_matrix(SNPmatrix & A, size_t c1, size_t c2, matrixType & M) {
  if(c1 >= A.size() || c2 >= A.size()) throw std::runtime_error("Bad bound in LD_matrix");
  const size_t n = c2-c1+1;
  if(n != M.nrow() || n != M.ncol()) {
    throw std::runtime_error("dimension mismatch in LD_matrix");
    return;
  }

  auto f = LD_pair_f<Algo, scalar_t>{};
#pragma omp parallel for
  for(size_t i1 = 0; i1 < n; i1++) {
    size_t x1 = c1+i1;
    for(size_t i2 = 0; i2 <= i1; i2++) {
      size_t x2 = c1+i2;
      M(i2, i1) = f(A, x1, x2);
    }
  } 

  // symetriser
#pragma omp parallel for
  for(size_t i1 = 0; i1 < n; i1++) {
    for(size_t i2 = 0; i2 < i1; i2++) {
       M(i1, i2) = M(i2, i1);
    }
  }
}

/***********************************************************************
 *
 * les fonctions pour calculer le LD de A[ c1:c2 ] avec A[ d1:d2 ]
 * plusieurs cas énumérés
 * la dernière fonction fait le dispatching
 * [code essentiellement repris de gaston. Pas eu le courage de faire plus 
 *  élégant]
 *
 ***********************************************************************/

/**** DON'T CALL THESE FUNCTIONS DIRECTLY, USE THE FINAL FUNCTION WHICH DOES THE DISPATCHING ****/

// Intervalles c1 c2 et d1 d2 disjoints [sauf possiblement un point sur la diagonale, mais pas de calculs en double]
template<LDalgorithm Algo, typename scalar_t = double, typename matrixType>
void LD_chunk_0(SNPmatrix & A, size_t c1, size_t c2, size_t d1, size_t d2, matrixType & M) {
  if(c2-c1+1 != M.nrow() || d2-d1+1 != M.ncol()) 
    throw std::runtime_error("dimension mismatch in LD_chunk_0");

  auto f = LD_pair_f<Algo, scalar_t>{};
#pragma omp parallel for
  for(size_t x2 = d1; x2 <= d2; x2++) {
    for(size_t x1 = c1; x1 <= c2; x1++) {
      M(x1 - c1, x2 - d1) = f(A, x1, x2);
    }
  }
}

// c1 <= d1 < c2 <= d2
template<LDalgorithm Algo, typename scalar_t = double, typename matrixType>
void LD_chunk_1(SNPmatrix & A, size_t c1, size_t c2, size_t d1, size_t d2, matrixType & M) {
  if(c2-c1+1 != M.nrow() || d2-d1+1 != M.ncol()) 
    throw std::runtime_error("dimension mismatch in LD_chunk_1");

  auto f = LD_pair_f<Algo, scalar_t>{};
#pragma omp parallel for
  for(size_t x2 = d1; x2 <= d2; x2++) 
    for(size_t x1 = c1; x1 < d1; x1++) 
      M(x1 - c1, x2 - d1) = f(A, x1, x2);

#pragma omp parallel for
  for(size_t x2 = d1; x2 <= c2; x2++) 
    for(size_t x1 = d1; x1 <= x2; x1++) 
      M(x1 - c1, x2 - d1) = f(A, x1, x2);

  // symetriser ce morceau
#pragma omp parallel for
  for(size_t x1 = d1; x1 <= c2; x1++) 
    for(size_t x2 = d1; x2 < x1; x2++) 
      M(x1 - c1, x2 - d1) = M(x2 - c1, x1 - d1);

#pragma omp parallel for
  for(size_t x2 = c2+1; x2 <= d2; x2++) 
    for(size_t x1 = d1; x1 <= c2; x1++) 
      M(x1 - c1, x2 - d1) = f(A, x1, x2);
}


// d1 <= c1 < d2 < =c2
template<LDalgorithm Algo, typename scalar_t = double, typename matrixType>
void LD_chunk_2(SNPmatrix & A, size_t c1, size_t c2, size_t d1, size_t d2, matrixType & M) {
  if(c2-c1+1 != M.nrow() || d2-d1+1 != M.ncol()) 
    throw std::runtime_error("dimension mismatch in LD_chunk_2");

  auto f = LD_pair_f<Algo, scalar_t>{};
#pragma omp parallel for
  for(size_t x2 = d1; x2 < c1; x2++) 
    for(size_t x1 = c1; x1 <= c2; x1++) 
      M(x1 - c1, x2 - d1) = f(A, x1, x2);
  
#pragma omp parallel for
  for(size_t x2 = c1; x2 <= d2; x2++) 
    for(size_t x1 = c1; x1 <= x2; x1++) 
      M(x1 - c1, x2 - d1) = f(A, x1, x2);
  
  // symetriser ce morceau
#pragma omp parallel for
  for(size_t x1 = c1; x1 <= d2; x1++) 
    for(size_t x2 = c1; x2 < x1; x2++) 
      M(x1 - c1, x2 - d1) = M(x2 - c1, x1 - d1);

#pragma omp parallel for
  for(size_t x2 = c1; x2 <= d2; x2++) 
    for(size_t x1 = d2+1; x1 <= c2; x1++) 
      M(x1 - c1, x2 - d1) = f(A, x1, x2);
}

// d1 < c1 <= c2 < d2
template<LDalgorithm Algo, typename scalar_t = double, typename matrixType>
void LD_chunk_3(SNPmatrix & A, size_t c1, size_t c2, size_t d1, size_t d2, matrixType & M) {
  if(c2-c1+1 != M.nrow() || d2-d1+1 != M.ncol()) 
    throw std::runtime_error("dimension mismatch in LD_chunk_3");

  auto f = LD_pair_f<Algo, scalar_t>{};
#pragma omp parallel for
  for(size_t x2 = d1; x2 < c1; x2++) 
    for(size_t x1 = c1; x1 <= c2; x1++) 
      M(x1 - c1, x2 - d1) = f(A, x1, x2);
  
#pragma omp parallel for
  for(size_t x2 = c1; x2 <= c2; x2++) 
    for(size_t x1 = c1; x1 <= x2; x1++) 
      M(x1 - c1, x2 - d1) = f(A, x1, x2);

  // symmetriser ce morceau
#pragma omp parallel for
  for(size_t x1 = c1; x1 <= c2; x1++) 
    for(size_t x2 = c1; x2 < x1; x2++) 
      M(x1 - c1, x2 - d1) = M(x2 - c1, x1 - d1);

#pragma omp parallel for
  for(size_t x2 = c2+1; x2 <= d2; x2++) 
    for(size_t x1 = c1; x1 <= c2; x1++) 
      M(x1 - c1, x2 - d1) = f(A, x1, x2);
}


// c1 <= d1 <= d2 <= c2 
template<LDalgorithm Algo, typename scalar_t = double, typename matrixType>
void LD_chunk_4(SNPmatrix & A, size_t c1, size_t c2, size_t d1, size_t d2, matrixType & M) {
  if(c2-c1+1 != M.nrow() || d2-d1+1 != M.ncol()) 
    throw std::runtime_error("dimension mismatch in LD_chunk_4");

  auto f = LD_pair_f<Algo, scalar_t>{};
#pragma omp parallel for
  for(size_t x2 = d1; x2 <= d2; x2++) 
    for(size_t x1 = c1; x1 < d1; x1++) 
      M(x1 - c1, x2 - d1) = f(A, x1, x2);
  
#pragma omp parallel for
  for(size_t x2 = d1; x2 <= d2; x2++) 
    for(size_t x1 = d1; x1 <= x2; x1++) 
      M(x1 - c1, x2 - d1) = f(A, x1, x2);

  // symetriser ce morceau
#pragma omp parallel for
  for(size_t x1 = d1; x1 <= d2; x1++) 
    for(size_t x2 = d1; x2 < x1; x2++) 
      M(x1 - c1, x2 - d1) = M(x2 - c1, x1 - d1);
  
#pragma omp parallel for
  for(size_t x2 = d1; x2 <= d2; x2++) 
    for(size_t x1 = d2+1; x1 <= c2; x1++) 
      M(x1 - c1, x2 - d1) = f(A, x1, x2);
}


// Cette fonction fait le choix de la bonne fonction
template<LDalgorithm Algo, typename scalar_t = double, typename matrixType>
void LD_chunk(SNPmatrix & A, size_t c1, size_t c2, size_t d1, size_t d2, matrixType & M) {
  if(c2 <= d1 || d2 <= c1)
    LD_chunk_0<Algo, scalar_t>(A, c1, c2, d1, d2, M);
  else if(c1 <= d1 && c2 <= d2)
    LD_chunk_1<Algo, scalar_t>(A, c1, c2, d1, d2, M);
  else if(d1 <= c1 && d2 <= c2)
    LD_chunk_2<Algo, scalar_t>(A, c1, c2, d1, d2, M);
  else if(d1 <= c1 && c2 <= d2)
    LD_chunk_3<Algo, scalar_t>(A, c1, c2, d1, d2, M);
  else if(c1 <= d1 && d2 <= c2)
    LD_chunk_4<Algo, scalar_t>(A, c1, c2, d1, d2, M);
  else 
    throw std::runtime_error("Uncatched case, this should not happen");
}

#endif


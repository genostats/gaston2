#include "SNPvector.h"
#include "SNPmatrix.h"

  // Question : getSNP ne pourrait pas renvoyer une référence à un snp ? A-t-on jamais besoin de récupérer le shared pointeur ?
  // Question : Est-ce que la parallélisation dans la fonction LD est efficace ? Je pense que dans gaston ça ne l'était pas, j'avais laissé tomber
  //            [ update : vu les résultats du microbenchmark la réponse est non, c'est 2 fois plus lent : j'ai supprimé la ligne ]
  //            On pourrait paralléliser sur la boucle la plus extérieure dans les fonctions LD diverses plutôt
  // NOTE : LD n'utilise pas le système de choix du centrage (plutôt une bonne idée a priori le centrage "p" n'est pas une bonne idée)

// calcule le LD d'une pair de SNP
template<typename scalar_t = double>
inline scalar_t LD_pair(SNPmatrix & M, size_t i1, size_t i2) {
  SNPvector & snp1 = *(M.getSNP(i1));
  // snp1.compute_stats();  // TODO a déplacer dans la fonction de calcul du LD, en ne recalculant pas si c'est déjà calculé. 
                         // Je laisse comme ça mais ça ralentit énormément, les stats sont recalculées plein de fois. Il faut revoir un peu la gestion de mu, sigma, etc.
                         // NOTE : je ne me souviens plus pourquoi dans le calcul de sigma le nbre de NAs intervient dans le calcul ?!

  SNPvector & snp2 = *(M.getSNP(i2));
  // snp2.compute_stats();

  return snp1.LD<scalar_t>(snp2);
}


// remplit une matrice carrée de LD des SNPS i avec c1 <= i <= c2 
// doit être ok avec n'importe quelle classe de matrice qui a des membres nrow() ncol() et l'affectation par M(i,j) = ...
// (column major mode matrix)
 
template<typename scalar_t = double, typename matrixType>
void LD_matrix(SNPmatrix & A, size_t c1, size_t c2, matrixType & M) {
  if(c1 >= A.size() || c2 >= A.size()) throw std::runtime_error("Bad bound in LD_matrix");
  const size_t n = c2-c1+1;
  if(n != M.nrow() || n != M.ncol()) {
    throw std::runtime_error("dimension mismatch in LD_matrix");
    return;
  }

  for(size_t i1 = 0; i1 < n; i1++) {
    size_t x1 = c1+i1;
    for(size_t i2 = 0; i2 <= i1; i2++) {
      size_t x2 = c1+i2;
      M(i2, i1) = LD_pair<scalar_t>(A, x1, x2);
    }
  } 

  // symetriser
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
 * [code repris de gaston. Pas eu le courage de réfléchir si on 
 * peut faire plus simple]
 *
 ***********************************************************************/

/**** DON'T CALL THESE FUNCTIONS DIRECTLY, USE THE FINAL FUNCTION WHICH DOES THE DISPATCHING ****/

// Intervalles c1 c2 et d1 d2 disjoints [sauf possiblement un point sur la diagonale, mais pas de calculs en double]
template<typename scalar_t = double, typename matrixType>
void LD_chunk_0(SNPmatrix & A, size_t c1, size_t c2, size_t d1, size_t d2, matrixType & M) {
  const int n = c2-c1+1;
  const int m = d2-d1+1;
  if(n != M.nrow() || m != M.ncol()) 
    throw std::runtime_error("dimension mismatch in LD_chunk_0");

  for(int x2 = d1; x2 <= d2; x2++) {
    for(int x1 = c1; x1 <= c2; x1++) {
      M(x1 - c1, x2 - d1) = LD_pair<scalar_t>(A, x1, x2);
    }
  }
}

// c1 <= d1 < c2 <= d2
template<typename scalar_t = double, typename matrixType>
void LD_chunk_1(SNPmatrix & A, size_t c1, size_t c2, size_t d1, size_t d2, matrixType & M) {
  if(c2-c1+1 != M.nrow() || d2-d1+1 != M.ncol()) 
    throw std::runtime_error("dimension mismatch in LD_chunk_1");

  for(int x2 = d1; x2 <= d2; x2++) 
    for(int x1 = c1; x1 < d1; x1++) 
      M(x1 - c1, x2 - d1) = LD_pair<scalar_t>(A, x1, x2);

  for(int x2 = d1; x2 <= c2; x2++) 
    for(int x1 = d1; x1 <= x2; x1++) 
      M(x1 - c1, x2 - d1) = LD_pair<scalar_t>(A, x1, x2);

  // symetriser ce morceau
  for(int x1 = d1; x1 <= c2; x1++) 
    for(int x2 = d1; x2 < x1; x2++) 
      M(x1 - c1, x2 - d1) = M(x2 - c1, x1 - d1);

  for(int x2 = c2+1; x2 <= d2; x2++) 
    for(int x1 = d1; x1 <= c2; x1++) 
      M(x1 - c1, x2 - d1) = LD_pair<scalar_t>(A, x1, x2);
}


// d1 <= c1 < d2 < =c2
template<typename scalar_t = double, typename matrixType>
void LD_chunk_2(SNPmatrix & A, size_t c1, size_t c2, size_t d1, size_t d2, matrixType & M) {
  if(c2-c1+1 != M.nrow() || d2-d1+1 != M.ncol()) 
    throw std::runtime_error("dimension mismatch in LD_chunk_2");

  for(int x2 = d1; x2 < c1; x2++) 
    for(int x1 = c1; x1 <= c2; x1++) 
      M(x1 - c1, x2 - d1) = LD_pair<scalar_t>(A, x1, x2);
  
  for(int x2 = c1; x2 <= d2; x2++) 
    for(int x1 = c1; x1 <= x2; x1++) 
      M(x1 - c1, x2 - d1) = LD_pair<scalar_t>(A, x1, x2);
  
  // symetriser ce morceau
  for(int x1 = c1; x1 <= d2; x1++) 
    for(int x2 = c1; x2 < x1; x2++) 
      M(x1 - c1, x2 - d1) = M(x2 - c1, x1 - d1);

  for(int x2 = c1; x2 <= d2; x2++) 
    for(int x1 = d2+1; x1 <= c2; x1++) 
      M(x1 - c1, x2 - d1) = LD_pair<scalar_t>(A, x1, x2);
}

// d1 < c1 <= c2 < d2
template<typename scalar_t = double, typename matrixType>
void LD_chunk_3(SNPmatrix & A, size_t c1, size_t c2, size_t d1, size_t d2, matrixType & M) {
  if(c2-c1+1 != M.nrow() || d2-d1+1 != M.ncol()) 
    throw std::runtime_error("dimension mismatch in LD_chunk_3");

  for(int x2 = d1; x2 < c1; x2++) 
    for(int x1 = c1; x1 <= c2; x1++) 
      M(x1 - c1, x2 - d1) = LD_pair<scalar_t>(A, x1, x2);
  
  for(int x2 = c1; x2 <= c2; x2++) 
    for(int x1 = c1; x1 <= x2; x1++) 
      M(x1 - c1, x2 - d1) = LD_pair<scalar_t>(A, x1, x2);

  // symmetriser ce morceau
  for(int x1 = c1; x1 <= c2; x1++) 
    for(int x2 = c1; x2 < x1; x2++) 
      M(x1 - c1, x2 - d1) = M(x2 - c1, x1 - d1);

  for(int x2 = c2+1; x2 <= d2; x2++) 
    for(int x1 = c1; x1 <= c2; x1++) 
      M(x1 - c1, x2 - d1) = LD_pair<scalar_t>(A, x1, x2);
}


// c1 <= d1 <= d2 <= c2 
template<typename scalar_t = double, typename matrixType>
void LD_chunk_4(SNPmatrix & A, size_t c1, size_t c2, size_t d1, size_t d2, matrixType & M) {
  if(c2-c1+1 != M.nrow() || d2-d1+1 != M.ncol()) 
    throw std::runtime_error("dimension mismatch in LD_chunk_4");

  for(int x2 = d1; x2 <= d2; x2++) 
    for(int x1 = c1; x1 < d1; x1++) 
      M(x1 - c1, x2 - d1) = LD_pair<scalar_t>(A, x1, x2);
  
  for(int x2 = d1; x2 <= d2; x2++) 
    for(int x1 = d1; x1 <= x2; x1++) 
      M(x1 - c1, x2 - d1) = LD_pair<scalar_t>(A, x1, x2);

  // symetriser ce morceau
  for(int x1 = d1; x1 <= d2; x1++) 
    for(int x2 = d1; x2 < x1; x2++) 
      M(x1 - c1, x2 - d1) = M(x2 - c1, x1 - d1);
  
  for(int x2 = d1; x2 <= d2; x2++) 
    for(int x1 = d2+1; x1 <= c2; x1++) 
      M(x1 - c1, x2 - d1) = LD_pair<scalar_t>(A, x1, x2);
}


// Cette fonction fait le choix de la bonne fonction
template<typename scalar_t = double, typename matrixType>
void LD_chunk(SNPmatrix & A, size_t c1, size_t c2, size_t d1, size_t d2, matrixType & M) {
  if(c2 <= d1 || d2 <= c1)
    LD_chunk_0<scalar_t>(A, c1, c2, d1, d2, M);
  else if(c1 <= d1 && c2 <= d2)
    LD_chunk_1<scalar_t>(A, c1, c2, d1, d2, M);
  else if(d1 <= c1 && d2 <= c2)
    LD_chunk_2<scalar_t>(A, c1, c2, d1, d2, M);
  else if(d1 <= c1 && c2 <= d2)
    LD_chunk_3<scalar_t>(A, c1, c2, d1, d2, M);
  else if(c1 <= d1 && d2 <= c2)
    LD_chunk_4<scalar_t>(A, c1, c2, d1, d2, M);
  else 
    throw std::runtime_error("Uncatched case, this should not happen");
}

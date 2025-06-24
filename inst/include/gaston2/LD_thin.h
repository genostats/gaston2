#include "SNPvector.h"
#include "SNPmatrix.h"
#include "Datastruct.h"
#include "Column.h"
#include "LD.h"
#include <stdexcept>
#include "debug.h"


#ifndef _LD_thinning_
#define _LD_thinning_


// cette fonction utilise l'estimateur des moments, en float. On n'a pas jugé
// utile de templater pour permettre autre chose.
//
// le type boolVector est un vecteur de booléens (opérateurs size() et [])
// le type FUN est soit une fonction soit une classe avec un operateur()
//
//  threshold  : seuil de thinning (sur r^2)
//
// max_dist_bp : distance en bp au delà de laquelle on ne considère le LD nul sans le calculer
// max_dist_cM : distance en cM au delà de laquelle on ne considère le LD nul sans le calculer
// UNE SEULE DE CES DEUX DISTANCES DOIT ETRE NON NULLE. Celle qui est nulle ne sera pas prise en compte.
//
// which_keep  : vecteur modifié en place. Doit être initialisé avec une longueur de M.nbSNPs() 
//              les SNPs pour lequels c'est déjà à false seront ignorés : initialiser avec des true
//              pour considérer tous les SNPs
//  keep_left : fonctions qui sera appelée par keep_left(i, j) avec i < j. Si elle renvoie 'true'
//              on garde le SNP i, sinon le j.
//              exemple : [](size_t i, size_t j) { return true; } pour toujours garder celui de gauche...
//
template<typename SNPvectorClass, typename boolVector, typename FUN>
void LD_thin(SNPmatrix<SNPvectorClass> & A, float threshold, unsigned int max_dist_bp, float max_dist_cM, boolVector which_keep, FUN keep_left) {

  size_t nbSNPs = A.nbSNPs();

  if(which_keep.size() != nbSNPs) 
    throw std::runtime_error("In ld_thin, dimensions of which_keep and A mismatch");

  // which distance will be used
  bool use_bp;
  if(max_dist_bp == 0) {
    if(max_dist_cM == 0) throw std::runtime_error("In ld_thin, both distances to 0");
    use_bp = false;
  } else {
    if(max_dist_cM != 0) throw std::runtime_error("In ld_thin, both distances non zero");
    use_bp = true;
  }

  const std::vector<int> & chr     = *A.getSNPStats().getColumn("chr").template get<int>();
  const std::vector<int> & pos     = *A.getSNPStats().getColumn("pos").template get<int>();
  const std::vector<double> & dist = *A.getSNPStats().getColumn("dist").template get<double>();
  // TODO test length of chr, pos, dist ??
  
  // we start at the beginning...
  size_t i = 0;

  threshold = std::sqrt(threshold); // on va calculer r, pas r^2
 
  while(i < nbSNPs) {
    size_t j = i + 1;
    int chr_i = chr[i];

    size_t next_i;
    bool found_next_i = false;

    int max_pos = pos[i] + max_dist_bp;
    double max_dist = dist[i] + max_dist_cM;

    // object which computed the LD
    auto ldf = LD_pair_f<LDalgorithm::moments, float, SNPvectorClass>{};

    while(j < nbSNPs && chr[i] == chr_i) {
      // test dist/pos according to use_bp...
      if(use_bp) {
        if(pos[j] > max_pos) break;
      } else {
        if(dist[j] > max_dist) break;
      }
      if(!which_keep[j]) { // SNP already removed
        j++;
        continue;
      }
      float ld_ij = ldf(A, i, j);
      if( std::abs(ld_ij) <= threshold ) { // LD is below threshold
        if(!found_next_i) {
          next_i = j;
          found_next_i = true;
        }
        j++;
      } else { // LD is above threshold
        if(keep_left(i,j)) { // keep i, remove j
          which_keep[j] = false;
          j++; // next loop will consider next j value
        } else { // remove i, keep j
          which_keep[i] = false;
          break; // break while(j) loop, will consider next i value
        }
      }
    } // end of while(j) loop
    if(found_next_i) i = next_i; else i = j;
  } // end of while(i) loop
}
#endif


#include "hwe.h"
#include "DataStruct.h"
#include "isAutosome.h"


#ifndef _gaston2_set_hwe_
#define _gaston2_set_hwe_

// prend un data struct "snp stats" et y ajoute une colonne hwe 
template<HWEtest test, typename scalar_t = double>
void set_hwe(DataStruct & snp_stats) {
  int nr = snp_stats.nrow();

  const std::vector<int> & N0 = *snp_stats.getColumn("N0").get<int>();
  const std::vector<int> & N1 = *snp_stats.getColumn("N1").get<int>();
  const std::vector<int> & N2 = *snp_stats.getColumn("N2").get<int>();

  std::vector<scalar_t> hwe;
  
  auto f = hwe_f<test, scalar_t>{};

  // y a-t-il des snp sur le chr X? si oui il y a les col N0.f N1.f N2.f
  if(snp_stats.hasColumn("N0.f")) {
    const std::vector<int> & chr = *snp_stats.getColumn("chr").get<int>();
    const std::vector<int> & N0f = *snp_stats.getColumn("N0.f").get<int>();
    const std::vector<int> & N1f = *snp_stats.getColumn("N1.f").get<int>();
    const std::vector<int> & N2f = *snp_stats.getColumn("N2.f").get<int>();
    for(int i = 0; i < nr; i++) {
      if(isX( chr[i] )) {
        hwe.push_back( f(N0f[i], N1f[i], N2f[i]) );
      } else {
        hwe.push_back( f(N0[i], N1[i], N2[i]) );
      }
    }
  } else {
    for(int i = 0; i < nr; i++) hwe.push_back( f(N0[i], N1[i], N2[i]) );
  }

  // ajout de la colonne
  snp_stats.setColumn(Column(hwe), "hwe");
}

#endif

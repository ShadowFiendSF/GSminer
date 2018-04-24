#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector elemLen1(const List& l)
{
  IntegerVector res;
  for(size_t i = 0; i != static_cast<size_t>(l.size()); ++i)
  {
    if(Rf_isNull(l[i]))
    {
      res.push_back(0L);
    }else{
      CharacterVector tmp = l[i];
      res.push_back(tmp.size());
    }
  }
  res.names() = l.names();
  return res;
}

#include "Rcpp.h"
#include <vector>
#include <cstdio>
#include <string>
#include <iostream>
using namespace Rcpp;
// [[Rcpp::export]]
RcppExport List GO2gene( const DataFrame &Sdata)
{
//	size_t GOsize = GOList.size();
	CharacterVector id = Sdata[1];
	CharacterVector uid = unique(id);
	size_t idSize = uid.size();
	if(idSize ==0)
	{
		return List::create();
	}
	List res(idSize);
	res.names() = uid;
	CharacterVector geneName(Sdata[0]);
	for(size_t i=0; i < static_cast<size_t>(idSize); ++i)
	{
		CharacterVector geneList;
//		std::string go(uid[i]);
		for(size_t j=0; j < static_cast<size_t>(id.size()); ++j)
		{
//			std::cout << uid[i] << "==" << geneName[j] <<std::endl;
			if(id[j] == uid[i]) geneList.push_back(geneName[j]); 
		}
		res[std::string(uid[i])] = unique(geneList);
	}
	return res;
}

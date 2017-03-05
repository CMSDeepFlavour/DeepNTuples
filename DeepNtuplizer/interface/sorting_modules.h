/*
 * sorting_modules.h
 *
 *  Created on: 5 Mar 2017
 *      Author: jkiesele
 */

#ifndef DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_SORTING_MODULES_H_
#define DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_SORTING_MODULES_H_

#include <algorithm>

namespace sorting{

template<class T>
bool comparePt(T a, T b){
	return a->pt()<b->pt();
}


template<class T>
bool compareDxy(T a, T b){
	return a->dxy()<b->dxy();
}

template<class T>
bool compareDxyDxyErr(T a, T b){
	float aerr=a->dxyError();
	float berr=b->dxyError();
	if( aerr && berr)
		return a->dxy()/aerr<b->dxy()/berr;
	if(!aerr || aerr!=aerr)
		return true;
	return false;
}


}
#endif /* DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_SORTING_MODULES_H_ */

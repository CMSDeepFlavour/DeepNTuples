/*
 * sorting_modules.h
 *
 *  Created on: 5 Mar 2017
 *      Author: jkiesele
 */

#ifndef DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_SORTING_MODULES_H_
#define DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_SORTING_MODULES_H_

#include <algorithm>
#include <cmath>
namespace sorting{

/*
 * the std::sort function needs a strict ordering.
 * means if your_function(a,b) is true, than your_function(b,a) must be false for all a and b
 */

template<class T>
bool comparePt(T a, T b){
	if(!a)return true;
	if(!b)return false;
	return a->pt()<b->pt();
}


template<class T>
bool compareDxy(T b, T a){
	if(!a && b)return true;
	if(!b && a)return false;
	if(!a && !b)return false;
	return a->dxy()<b->dxy();
}

template<class T>
bool compareDxyDxyErr(T b, T a){
	if(!a) return true;
	if(!b) return false;
	if(!a && !b)return false;

	float aerr=a->dxyError();
	float berr=b->dxyError();


	float asig=a->dxy()/aerr;
	float bsig=b->dxy()/berr;

	if(!std::isnormal(asig) && std::isnormal(bsig))
		return true;
	else if(!std::isnormal(bsig) && std::isnormal(asig))
		return false;
	else if(!std::isnormal(bsig) && !std::isnormal(asig))
		return false;

	return asig<bsig;
}


}
#endif /* DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_SORTING_MODULES_H_ */

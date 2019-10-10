/*
 * BaseParticle.cpp
 *
 *  Created on: 21/set/2010
 *      Author: lorenzo
 */

#include "BaseParticle.h"

#include "../Boxes/BaseBox.h"

template<typename number>
BaseParticle<number>::BaseParticle() : index(-1), type(P_INVALID), n3(P_VIRTUAL), n5(P_VIRTUAL) {
	en3 = (number) 0;
	en5 = (number) 0;
	esn3 = (number) 0;
	esn5 = (number) 0;
	inclust = false;
	ext_potential = (number) 0.;
	strand_id = -1;
	N_int_centers = 0;
	//_pos_shift = LR_vector<number>(0., 0., 0.);
	_pos_shift[0] = 0; _pos_shift[1] = 0; _pos_shift[2] = 0;
	force = LR_vector<number>(0., 0., 0.);
	torque = LR_vector<number>(0., 0., 0.);
	int_centers = NULL;
	btype = 0;
	next_particle = P_INVALID;
}

template<typename number>
void BaseParticle<number>::copy_from(const BaseParticle<number> &p) {
	index = p.index;
	type = p.type;
	btype = p.btype;
	pos = p.pos;
	vel = p.vel;
	orientation = p.orientation;
	orientationT = p.orientationT;
	pos_list = p.pos_list;
	force = p.force;
	en3 = p.en3;
	en5 = p.en5;
	esn3 = p.esn3;
	esn5 = p.esn5;
	n3 = p.n3;
	n5 = p.n5;

	for(int i = 0; i < N_int_centers; i++) int_centers[i] = p.int_centers[i];

	ext_potential = p.ext_potential;
}

template<typename number>
BaseParticle<number>::~BaseParticle() {
	if(int_centers != NULL) delete[] int_centers;
}

template<typename number>
void BaseParticle<number>::set_initial_forces (llint step, const BoxPtr<number> &box) {
	LR_vector<number> abs_pos = box->get_abs_pos(this);
	if (this->is_rigid_body()) this->torque = LR_vector<number>((number)0.f, (number)0.f, (number)0.f);
	this->force = LR_vector<number>((number)0.f, (number)0.f, (number)0.f);
	for(auto ext_force : ext_forces) {
		this->force += ext_force->value(step, abs_pos);
	}
}

template<typename number>
bool BaseParticle<number>::add_ext_force(ForcePtr<number> f) {
	ext_forces.push_back(f);

	return true;
}

template<typename number>
void BaseParticle<number>::init() {
	force = LR_vector<number>(0., 0., 0.);
	torque = LR_vector<number>(0., 0., 0.);
	_check();
}

template<typename number>
void BaseParticle<number>::_check() {
	assert(index >= 0);
	assert(type != P_INVALID);
}

template class BaseParticle<double>;
template class BaseParticle<float>;

#include "Clock.h"

Clock::Clock(double start_time, double end_time, double dt_outer, double dt_inner){

  /* Initialize variables */
  _times[START] = 0.0;
  _times[PREVIOUS_OUT] = 0.0;
  _times[PREVIOUS_IN] = 0.0;
  _times[CURRENT] = 0.0;
  _times[FORWARD_OUT] = 0.0;
  _times[END] = 0.0;

  /* Set time step sizes */
  setDtInner(dt_inner);
  setDtOuter(dt_outer);
  setStartTime(start_time);
  setEndTime(end_time);
}


Clock::~Clock(){

  _times.clear();
}


double Clock::getTime(int state){

  if (state == FORWARD_IN_OLD)
    return _times[CURRENT];
  else if (state == FORWARD_OUT_OLD)
    return _times[FORWARD_OUT];
  else
    return _times[state];
}


double Clock::getDtInner(){
  return _dt_inner;
}


double Clock::getDtOuter(){
  return _dt_outer;
}


void Clock::setDtInner(double dt_inner){
  _dt_inner = dt_inner;
}


void Clock::setDtOuter(double dt_outer){
  _dt_outer = dt_outer;
}


void Clock::setTime(int state, double time){
  _times[state] = time;
}


void Clock::setStartTime(double time){
  _times[START] = time;
}


void Clock::setEndTime(double time){
  _times[END] = time;
}


void Clock::takeInnerStep(){
  _times[PREVIOUS_IN] = _times[CURRENT];
  _times[CURRENT] += _dt_inner;
}


void Clock::takeOuterStep(){

  _times[PREVIOUS_OUT] = _times[FORWARD_OUT];
  _times[FORWARD_OUT] = _times[FORWARD_OUT] + _dt_outer;  
  
  if (_times[FORWARD_OUT] > _times[END])
    _times[FORWARD_OUT] = _times[END];

  _times[CURRENT] = _times[PREVIOUS_OUT];
  _times[PREVIOUS_IN] = _times[PREVIOUS_OUT];
}


void Clock::resetToPreviousOuterStep(){
  _times[CURRENT] = _times[PREVIOUS_OUT];
  _times[PREVIOUS_IN] = _times[PREVIOUS_OUT];
}


std::string Clock::toString(){

  std::stringstream string;
  string << std::setprecision(6);
    
  string << "Clock" << std::endl;
  string << " START\t\t\t = " << _times[START] << std::endl;
  string << " PREVIOUS_OUT\t = " << _times[PREVIOUS_OUT] << std::endl;
  string << " PREVIOUS_IN\t = " << _times[PREVIOUS_IN] << std::endl;
  string << " CURRENT\t\t = " << _times[CURRENT] << std::endl;
  string << " FORWARD_OUT\t = " << _times[FORWARD_OUT] << std::endl;
  string << " END\t\t\t\t = " << _times[END] << std::endl;

  return string.str();
}


void Clock::printString(){
  log_printf(NORMAL, toString().c_str());
}


std::string Clock::getStateName(int state){

  std::string name;

  if (state == START)
    name = "START";
  else if (state == PREVIOUS_OUT)
    name = "PREVIOUS_OUT";
  else if (state == PREVIOUS_IN)
    name = "PREVIOUS_IN";
  else if (state == CURRENT)
    name = "CURRENT";
  else if (state == FORWARD_OUT)
    name = "FORWARD_OUT";
  else if (state == FORWARD_OUT_OLD)
    name = "FORWARD_OUT_OLD";
  else if (state == FORWARD_IN_OLD)
    name = "FORWARD_IN_OLD";
  else if (state == END)
    name = "END";
  else
    log_printf(ERROR, "Unable to get clock state name with state %d",
               state);

  return name;
}

#include "Clock.h"

Clock::Clock(double start_time, double end_time, double dt_outer, double dt_inner){

  /* Initialize variables */
  _times[START] = 0.0;
  _times[PREVIOUS_OUT] = 0.0;
  _times[PREVIOUS_IN] = 0.0;
  _times[CURRENT] = 0.0;
  _times[FORWARD_IN_OLD] = 0.0;
  _times[FORWARD_OUT] = 0.0;
  _times[FORWARD_OUT_OLD] = 0.0;
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


double Clock::getTime(int position){
  return _times[position];
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


void Clock::setTime(int position, double time){
  _times[position] = time;
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
  string << " FORWARD_IN_OLD\t = " << _times[FORWARD_IN_OLD] << std::endl;
  string << " FORWARD_OUT\t = " << _times[FORWARD_OUT] << std::endl;
  string << " FORWARD_OUT_OLD\t = " << _times[FORWARD_OUT_OLD] << std::endl;  
  string << " END\t\t\t\t = " << _times[END] << std::endl;

  return string.str();
}


void Clock::printString(){
  log_printf(NORMAL, toString().c_str());
}


std::string Clock::getPositionName(int position){

  std::string name;

  if (position == START)
    name = "START";
  else if (position == PREVIOUS_OUT)
    name = "PREVIOUS_OUT";
  else if (position == PREVIOUS_IN)
    name = "PREVIOUS_IN";
  else if (position == CURRENT)
    name = "CURRENT";
  else if (position == FORWARD_IN_OLD)
    name = "FORWARD_IN_OLD";
  else if (position == FORWARD_OUT)
    name = "FORWARD_OUT";
  else if (position == FORWARD_OUT_OLD)
    name = "FORWARD_OUT_OLD";
  else if (position == END)
    name = "END";
  else
    log_printf(ERROR, "Unable to get clock position name with position %d",
               position);

  return name;
}

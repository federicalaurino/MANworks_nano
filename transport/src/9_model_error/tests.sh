
echo ///////////////A REF/////////////////// 2>&1 | tee -a buffer_scenario.hpp

./M3D1D input_ScenarioA_REF.param 2>&1 | tee -a buffer_scenario.hpp

echo ///////////////A RED/////////////////// 2>&1 | tee -a buffer_scenario.hpp

./M3D1D input_ScenarioA_RED.param 2>&1 | tee -a buffer_scenario.hpp

echo ///////////////B REF/////////////////// 2>&1 | tee -a buffer_scenario.hpp

./M3D1D input_ScenarioB_REF.param 2>&1 | tee -a buffer_scenario.hpp

echo ///////////////B RED/////////////////// 2>&1 | tee -a buffer_scenario.hpp

./M3D1D input_ScenarioB_RED.param 2>&1 | tee -a buffer_scenario.hpp

echo ///////////////C REF/////////////////// 2>&1 | tee -a buffer_scenario.hpp

./M3D1D input_ScenarioC_REF.param 2>&1 | tee -a buffer_scenario.hpp

echo ///////////////C RED///////////////////2>&1 | tee -a buffer_scenario.hpp

./M3D1D input_ScenarioC_RED.param 2>&1 | tee -a buffer_scenario.hpp

echo ///////////////D REF/////////////////// 2>&1 | tee -a buffer_scenario.hpp

./M3D1D input_ScenarioD_REF.param 2>&1 | tee -a buffer_scenario.hpp

echo ///////////////D RED/////////////////// 2>&1 | tee -a buffer_scenario.hpp

./M3D1D input_ScenarioD_RED.param 2>&1 | tee -a buffer_scenario.hpp

echo ///////////////E REF/////////////////// 2>&1 | tee -a buffer_scenario.hpp

./M3D1D input_ScenarioE_REF.param 2>&1 | tee -a buffer_scenario.hpp
#./M3D1D input_ScenarioC_RED_uniform.param 2>&1 | tee -a buffer_scenario.hpp
#./M3D1D input_ScenarioC_REF_uniform.param 2>&1 | tee -a buffer_scenario.hpp


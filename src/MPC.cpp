#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

// TODO: Set the timestep length and duration
size_t N = 15;
double dt = 0.08;
// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;
const int COST_IDX = 0;
const int X_START = 0;
const int Y_START = X_START + N;
const int PSI_START = Y_START + N;
const int V_START = PSI_START + N;
const int CTE_START = V_START + N;
const int EPSI_START = CTE_START + N;
const int DELTA_START = EPSI_START + N;
const int ACC_START = DELTA_START + N - 1;
const double REF_V = 15.0;



class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    // TODO: implement MPC
    // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
    // NOTE: You'll probably go back and forth between this function and
    // the Solver function below.
    // DEFINE COST FUNCTION AND CONSTRAIN
    fg[0] = 0;
    // COST FUNCTION
    // add errors from CTE and steady speed 
    for (int i = 0; i < N; i++){
        fg[0] +=  5000 * CppAD::pow(vars[CTE_START + i], 2);
        fg[0] +=  5000 * CppAD::pow(vars[EPSI_START+ i], 2);
        fg[0] += 2 * CppAD::pow(vars[V_START+ i] - REF_V, 2); 
        //fg[0] += 20000 * CppAD::pow(3 * coeffs[2] * CppAD::pow(vars[X_START +i], 2), 2);
 
    }
    
    // actuators need to be smooth so we are not acc fast off the path 
    for (int i = 0; i < N - 1; i++){
        fg[0] += 120000 * CppAD::pow(vars[DELTA_START + i], 2);
        fg[0] += 30 * CppAD::pow(vars[ACC_START+ i], 2);
        // make sure we dont turn too steep while speed up 
        fg[0] += 100000 * CppAD::pow(vars[V_START+ i] - vars[V_START+ i + 1], 2) * vars[DELTA_START + i];
    }

    // minimize gap between actuation state so changing lane will be smoother
    for (int i = 0; i < N - 2; i++){
        fg[0] += 1000000 * CppAD::pow(vars[DELTA_START + i] - vars[DELTA_START + i + 1], 2);
        fg[0] += CppAD::pow(vars[ACC_START+ i] - vars[ACC_START+ i + 1], 2);
         
    }

    
    // INITIALIZE CONSTRAINS
    fg[X_START + 1] = vars[X_START];
    fg[Y_START + 1] = vars[Y_START];
    fg[PSI_START + 1] = vars[PSI_START];
    fg[V_START + 1] = vars[V_START];
    fg[CTE_START + 1] = vars[CTE_START];
    fg[EPSI_START + 1] = vars[EPSI_START];
    
   
    // set on subsequent constrain based on model state update equations
    for (int i = 1; i < N; i++){
         AD<double> x0 = vars[X_START + i -1];
         AD<double> f0 = coeffs[0] * coeffs[1] * x0 + coeffs[2] * x0 * x0 + coeffs[3] * x0 * x0 * x0;
         // expected epsi from model
         AD<double> epsi_pred = CppAD::atan(3 * coeffs[3] * x0 * x0 + 2 * coeffs[2] * x0 + coeffs[1]);        

         fg[X_START + i + 1] = vars[X_START + i] - 
            (vars[X_START + i - 1] + vars[V_START + i - 1]* CppAD::cos(vars[PSI_START + i - 1]) * dt);
        
        fg[Y_START + i + 1] = vars[Y_START + i] - 
            (vars[Y_START + i - 1] + vars[V_START + i - 1] * CppAD::sin(vars[PSI_START + i - 1]) * dt);
        
        fg[PSI_START + i + 1] = vars[PSI_START + i] -
                    (vars[PSI_START + i - 1] + vars[V_START + i - 1] * vars[DELTA_START + i - 1] * dt / Lf );  
        
        fg[V_START + i + 1] = vars[V_START + i] - (vars[V_START + i - 1] + vars[ACC_START + i - 1] * dt);

        fg[CTE_START + i + 1] = vars[CTE_START + i] - 
                (f0 - vars[Y_START + i - 1] + vars[V_START + i - 1] * CppAD::sin(vars[EPSI_START + i - 1]) * dt);
        
        fg[EPSI_START + i + 1] = vars[EPSI_START + i] -
                     (vars[PSI_START + i - 1] - epsi_pred 
                      + vars[V_START + i - 1] * vars[DELTA_START + i - 1] / Lf * dt);

    }




  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;
  size_t i;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  // TODO: Set the number of model variables (includes both states and inputs).
  // For example: If the state is a 4 element vector, the actuators is a 2
  // element vector and there are 10 timesteps. The number of variables is:
  //
  // 4 * 10 + 2 * 9
  size_t n_vars = 6 * N + 2 * (N -1);
  // TODO: Set the number of constraints
  size_t n_constraints = N * 6;

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }
  
  // initialize first value of variable
   vars[X_START] = state[0];
   vars[Y_START] = state[1];
   vars[PSI_START] = state[2];
   vars[V_START] = state[3];
   vars[CTE_START] = state[4];
   vars[EPSI_START] = state[5];



  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  // TODO: Set lower and upper limits for variables.
  for (int i = X_START; i < DELTA_START; i++) {
      vars_lowerbound[i] = -1.0e19;  
      vars_upperbound[i] = 1.0e19;
  }
 
  // Steering limitations
   for (int i = DELTA_START; i < ACC_START; i++) {
        vars_lowerbound[i] = -0.436332;
        vars_upperbound[i] = 0.436332;
  }

  // acceleration limitations
  for (int i = ACC_START; i < n_vars; i++) {
        vars_lowerbound[i] = -1.0;
        vars_upperbound[i] = 1.0;
  }


  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (int i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }

  // set constrain for initial value 
  constraints_lowerbound[X_START] = state[0];
  constraints_lowerbound[Y_START] = state[1];
  constraints_lowerbound[PSI_START] = state[2];
  constraints_lowerbound[V_START] = state[3];
  constraints_lowerbound[CTE_START] = state[4];
  constraints_lowerbound[EPSI_START] = state[5];

  constraints_upperbound[X_START] = state[0];
  constraints_upperbound[Y_START] = state[1];
  constraints_upperbound[PSI_START] = state[2];
  constraints_upperbound[V_START] = state[3];
  constraints_upperbound[CTE_START] = state[4];
  constraints_upperbound[EPSI_START] = state[5];



  // object that computes objective and constraints
  FG_eval fg_eval(coeffs);

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;
  std::cout << "Solution Size : " << solution.x.size() << std::endl;

  // TODO: Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  //
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.
  vector<double> result = {solution.x[DELTA_START], solution.x[ACC_START]};
  
  for (int i = X_START; i < Y_START; i++){
     result.push_back(solution.x[i]);
  }

  for (int i = Y_START; i < PSI_START; i++){
     result.push_back(solution.x[i]);
  }

  return result; 
}

## Model Predictive Control (MPC) Project Writeup

### Model Details
The MPC contains of 6 states:
    - cartesian coordinates x, y which respects to the orientation of the machine.
    - heading angle (psi) with respect to the x coordinate.
    - velocity of the car
    - cross track error which is the distance from the car coordinate to its planned trajectory.
    - orientation error which is the error of the car heading angle compare with its planned heading based on trajectory.

The MPC MPC also contains 2 control commands:
    - steering angle
    - acceleration 

The actuators will be optimized using a non-linear optimizer IPOPT with constrained on the below update equations:
<image src="model.png"/> 



### Timestep Lenght and Ealpsed Duration
Starting with N = 10 and dt = 0.1, the car only run fine for a small portion of the route but having a difficult time to drive close to the planned trajectory.
In order to help the ipopt solver find a better fit trajectory, after some trials and errors, increaing the number of timestep length N to 20 as well as lowering dt 0.08 to accomodate for actuated latency (100ms) seems to make the car to stabilize better and trajectory look smoother.

### Polynomial Fitting and MPC Perprocessing
All the waypoints recevied from simulator are with respected to global map coordinates. 
Thus they are needed to be transform to the car coordinates using translation and transformation matrix as below:

         for (int i = 0; i < ptsx.size(); i++) {
             double translation_x = ptsx[i] - px;
             double translation_y = ptsy[i] - py;
             // apply rotation here 
             //std::cout << "psi " << psi << " trans " << translation_x << std::endl;
             way_x.push_back(translation_x * cos(psi) + (translation_y) * sin(psi));
             way_y.push_back((-1) * (translation_x * sin(psi)) + translation_y * cos(psi));

          }
 

Also, we need to calculate the initial orientation erorr and crosstrack error with respect to car coordinate system as below:

Finally, we feed the initial state to the solver to optimize the actuators command.


### Model Predictive Control with Latency
In order to due with control latency, I first reduce the dt (elapsed between timesteps) to 0.08 so we can accomodate for some latency. Also, after observing the car is not slowing down soon enough to make the tight left turn after crossing the bridge, I add another penalty term which related to the product of veloctiy and steering angle at a specific moment in order to help the car slow down at steep turn and not running off the road due to some latency.




